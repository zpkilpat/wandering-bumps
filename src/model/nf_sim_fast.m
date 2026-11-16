function out = nf_sim_fast(varargin)
%NF_SIM_FAST  Reduced-order solver for the stochastic ring field (1.7), no input.
%
%   With no input (I0 = 0) the weight and the noise only drive the cos x and
%   sin x modes, so the field is exactly
%        U(x,t) = a(t) cos x + b(t) sin x = R cos(x - phi),
%   and the active set {U >= theta} is one arc |x - phi| <= alpha,
%   alpha = acos(theta/R) (empty if R < theta).  The integrals in (1.7),
%        c1 = int cos y (1 + sigma cos ny) H dy,   c2 = int sin y (1 + sigma cos ny) H dy,
%   then have closed forms, and each realization is the 2D SDE
%        da = (-a + c1) dt + sqrt(pi eps dt) xi1,   db = (-b + c2) dt + sqrt(pi eps dt) xi2.
%   Same Euler-Maruyama step and time step as NF_SIM, but no spatial grid, so
%   it solves the continuum model.  It is 50-100x faster than the grid solver
%   and agrees with it statistically (run  fast_vs_grid_check  to see).
%
%   Same name-value options and output struct as NF_SIM; I0 must be 0 (use the
%   grid solver for input-driven runs).  The bump centre is phi (the peak).
    p = struct('theta',0.5,'eps',0.01,'n',1,'sigma',0,'I0',0,'nI',1, ...
               'T',50,'dt',0.01,'M',1,'N',628,'rec',10,'seed',[],'relax',30, ...
               'keepfield',false,'stop','none','stopthr',3*pi/4,'fast',true);
    for k = 1:2:numel(varargin), p.(varargin{k}) = varargin{k+1}; end
    if p.I0 ~= 0, error('nf_sim_fast:input','Fast mode needs I0 = 0; use nf_sim.'); end
    if ~isempty(p.seed), rng(p.seed); end
    th = p.theta;  n = p.n;  sg = p.sigma;
    N = p.N;  dx = 2*pi/N;  x = -pi + (0:N-1)'*dx;
    sq = sqrt(p.eps*pi*p.dt);

    % stable initial bump, relaxed deterministically
    a = sqrt(1+th) + sqrt(1-th);  b = 0;
    for k = 1:round(p.relax/p.dt)
        [c1,c2] = drive(a,b,th,n,sg);
        a = a + p.dt*(-a + c1);  b = b + p.dt*(-b + c2);
    end
    c0 = atan2(b,a);
    a0 = a;  b0 = b;

    M = p.M;  nsteps = round(p.T/p.dt);  nrec = floor(nsteps/p.rec) + 1;
    a = repmat(a0,1,M);  b = repmat(b0,1,M);
    idx = 1:M;
    X = NaN(M,nrec);  event = NaN(M,1);
    field = [];  if p.keepfield, field = NaN(N,nrec); end
    t = (0:nrec-1)*p.rec*p.dt;

    for k = 0:nsteps
        if mod(k,p.rec) == 0
            r = k/p.rec + 1;
            ext = hypot(a,b) < th;
            rel = angle(exp(1i*(atan2(b,a) - c0)));
            X(idx,r) = rel(:);
            if r > 1                                        % extinct bump: hold last position
                X(idx(ext),r) = X(idx(ext),r-1);
            end
            if p.keepfield && idx(1) == 1, field(:,r) = a(1)*cos(x) + b(1)*sin(x); end
            switch p.stop
                case 'extinct', done = ext;
                case 'switch',  error('nf_sim_fast:switch','Switch mode needs input; use nf_sim.');
                otherwise
                    done = false(1,numel(idx));
                    new = ext & isnan(event(idx).');
                    event(idx(new)) = t(r);
            end
            if any(done)
                event(idx(done)) = t(r);
                a = a(~done);  b = b(~done);  idx = idx(~done);
                if isempty(idx), break; end
            end
        end
        if k == nsteps, break; end
        m = numel(idx);
        [c1,c2] = drive(a,b,th,n,sg);
        a = a + p.dt*(-a + c1) + sq*randn(1,m);
        b = b + p.dt*(-b + c2) + sq*randn(1,m);
    end
    if strcmp(p.stop,'none'), X = unwrap(X,[],2); end
    out = struct('t',t,'X',X,'event',event,'field',field,'x',x,'c0',c0,'params',p);
end

function [c1,c2] = drive(a,b,th,n,sg)
% closed-form integrals of cos y (1+sg cos ny) and sin y (1+sg cos ny) over the arc
    R   = hypot(a,b);
    al  = acos(min(th./R,1));                 % half-width; 0 if R < th (empty set)
    ph  = atan2(b,a);
    sa  = sin(al);
    S1  = sin((n+1)*al)/(n+1);                % int-type factor for mode n+1
    if n == 1, S0 = al; else, S0 = sin((n-1)*al)/(n-1); end
    c1 = 2*cos(ph).*sa + sg*(cos((n+1)*ph).*S1 + cos((n-1)*ph).*S0);
    c2 = 2*sin(ph).*sa + sg*(sin((n+1)*ph).*S1 - sin((n-1)*ph).*S0);
end
