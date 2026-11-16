function out = nf_sim(varargin)
%NF_SIM  Euler-Maruyama simulation of the stochastic ring neural field (1.7)
%
%   dU(x,t) = [ -U + int_{-pi}^{pi} w(x,y) f(U(y,t)) dy + I(x) ] dt + sqrt(eps) dW(x,t)
%
%   f(u)   = H(u - theta)                                  Heaviside        (1.6)
%   w(x,y) = (1 + sigma cos(n y)) cos(x - y)               weight           (1.4)
%   I(x)   = I0 cos(nI x)                                  input            (2.3)
%   <dW(x,t) dW(y,t)> = C(x-y) dt,  C(x) = pi cos(x)       noise            (4.1)
%
%   The weight and the noise are both rank two in (cos, sin), so the
%   convolution and the noise are computed exactly in O(N) per realization:
%       int w f dy = cos(x) * int cos(y)(1+sigma cos ny) f dy
%                  + sin(x) * int sin(y)(1+sigma cos ny) f dy
%       dW         = sqrt(pi) [cos(x) dB1 + sin(x) dB2]
%   Integrals use the trapezoidal rule on the periodic grid (a plain sum * dx).
%
%   Name-value options (defaults follow the paper):
%     'theta' 0.5   'eps' 0.01   'n' 1   'sigma' 0   'I0' 0   'nI' 1
%     'T' 50   'dt' 0.01   'M' 1  (realizations)   'N' 628 (dx ~ 0.01)
%     'rec' 10 (record every rec steps)   'seed' []   'relax' 30
%     'keepfield' false   store U(x,t) of realization 1
%     'fast' []  true: use the exact two-variable reduction NF_SIM_FAST (about 40x faster;
%            only without input).  Default follows run_all_figures('fast',true).
%     'stop' 'none' | 'extinct' | 'switch'   stop and drop a realization at its
%            event; 'switch' means |center| first exceeds 'stopthr' (3 pi/4)
%
%   The initial condition is the stable bump: A cos x with the homogeneous
%   amplitude (A.8), relaxed deterministically for 'relax' time units so that
%   it sits exactly on the attractor of the heterogeneous / input-driven
%   network.  Positions are measured relative to that centre.
%
%   out.t        1 x nrec  recording times
%   out.X        M x nrec  bump centre relative to the start, unwrapped in
%                          time ('none' mode) or wrapped to (-pi,pi] ('extinct'/'switch')
%   out.event    M x 1     event time: extinction ('none','extinct') or switch
%                          ('switch'); NaN if it never happened
%   out.field    N x nrec  U(x,t) of realization 1 (if keepfield)
%   out.x        N x 1     grid
%   out.c0       start position of the relaxed bump
    p = struct('theta',0.5,'eps',0.01,'n',1,'sigma',0,'I0',0,'nI',1, ...
               'T',50,'dt',0.01,'M',1,'N',628,'rec',10,'seed',[],'relax',30, ...
               'keepfield',false,'stop','none','stopthr',3*pi/4,'fast',[]);
    for k = 1:2:numel(varargin), p.(varargin{k}) = varargin{k+1}; end
    % fast mode: exact 2D reduction (see NF_SIM_FAST); needs no input and no switch criterion
    if isempty(p.fast), p.fast = isappdata(0,'wb_fast') && getappdata(0,'wb_fast'); end
    if p.fast && p.I0 == 0 && ~strcmp(p.stop,'switch')
        out = nf_sim_fast(varargin{:});  return
    end
    if ~isempty(p.seed), rng(p.seed); end

    N = p.N;  dx = 2*pi/N;  x = -pi + (0:N-1)'*dx;
    cx = cos(x);  sx = sin(x);
    hx = 1 + p.sigma*cos(p.n*x);
    Ck = [cx.*hx, sx.*hx]' * dx;                    % 2 x N, includes trapezoid weight
    Inp = p.I0*cos(p.nI*x);
    sq = sqrt(p.eps*pi*p.dt);

    % --- stable initial bump, relaxed deterministically ------------------
    A  = sqrt(1+p.theta) + sqrt(1-p.theta);
    u0 = A*cx + Inp;
    for k = 1:round(p.relax/p.dt)
        f  = double(u0 >= p.theta);
        c  = Ck*f;
        u0 = u0 + p.dt*(-u0 + cx*c(1) + sx*c(2) + Inp);
    end
    f0 = double(u0 >= p.theta);
    c0 = atan2(sx'*f0, cx'*f0);

    % --- stochastic run --------------------------------------------------
    M = p.M;  nsteps = round(p.T/p.dt);  nrec = floor(nsteps/p.rec) + 1;
    u = repmat(u0, 1, M);
    idx = 1:M;                                      % realizations still running
    X = NaN(M, nrec);  event = NaN(M,1);
    field = [];  if p.keepfield, field = NaN(N, nrec); end
    t = (0:nrec-1)*p.rec*p.dt;
    stopmode = p.stop;

    for k = 0:nsteps
        if mod(k,p.rec) == 0
            r = k/p.rec + 1;
            f = double(u >= p.theta);
            ang = atan2(sx'*f, cx'*f);                  % centre of the active set
            rel = angle(exp(1i*(ang - c0)));
            X(idx, r) = rel(:);
            if r > 1                                    % an extinct bump has no centre:
                dead = (sum(f,1) == 0);                 % hold its last position
                X(idx(dead), r) = X(idx(dead), r-1);
            end
            if p.keepfield && idx(1) == 1, field(:,r) = u(:,1); end
            ext = (sum(f,1) == 0);                       % bump extinct (no active cells)
            switch stopmode
                case 'extinct'
                    done = ext;
                case 'switch'
                    done = abs(rel) >= p.stopthr;
                otherwise
                    done = false(1,numel(idx));
                    new = ext & isnan(event(idx).');
                    event(idx(new)) = t(r);              % record extinction, keep running
            end
            if any(done)
                event(idx(done)) = t(r);
                u = u(:,~done);  idx = idx(~done);
                if isempty(idx), break; end
            end
        end
        if k == nsteps, break; end
        m = numel(idx);
        f = double(u >= p.theta);
        c = Ck*f;                                       % 2 x m
        u = u + p.dt*(-u + cx*c(1,:) + sx*c(2,:) + Inp) ...
              + sq*(cx*randn(1,m) + sx*randn(1,m));
    end

    if strcmp(stopmode,'none')
        X = unwrap(X, [], 2);                           % continuous position
    end
    out = struct('t',t,'X',X,'event',event,'field',field,'x',x,'c0',c0, ...
                 'params',p);
end
