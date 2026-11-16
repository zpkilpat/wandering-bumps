function th = wb_theory()
%WB_THEORY  Analytical formulas from Kilpatrick & Ermentrout, SIADS 12(1), 2013.
%
%   th = wb_theory() returns a struct of function handles.  Equation numbers
%   refer to the paper.  Ring model on [-pi,pi], Heaviside firing rate,
%   cosine kernel wbar(x) = cos(x), weight (1+sigma cos(n y)) cos(x-y),
%   noise correlation C(x) = pi cos(x).
%
%   Homogeneous network (Appendix A, sec. 4.1)
%     th.Ahom(theta,s)          bump amplitude A = sqrt(1+theta) + s sqrt(1-theta)   (A.8)
%     th.ahom(theta,s)          half-width  a = acos(theta/A)                        (A.9)
%     th.Dhom(eps,theta)        diffusion coefficient D(eps)                         (4.3)
%
%   Input-locked bumps (sec. 2.2, 4.3)   I(x) = I0 cos(n x)
%     th.a_input(I0,n,theta)    half-widths (narrow, wide) solving (2.6)
%     th.kappa_input(I0,n,a)    mean-reversion rate                                  (4.4)
%     th.ou_var(D,kappa,t)      Ornstein-Uhlenbeck variance                          (3.22)
%
%   Synaptic heterogeneity (sec. 2.3, 3.3, 3.4)   w = (1+sigma cos(n y)) cos(x-y)
%     th.Ahet(a,sg,n,sigma)     bump amplitude A_+- as a function of a               (2.15)
%     th.a_het(theta,sg,n,sigma) half-widths (narrow, wide)                          (2.16),(2.17)
%     th.lambda_odd(theta,sg,n,sigma)  odd (shift) eigenvalue, direct from           (2.18)
%     th.lambda_odd_closed(...)  closed form for n>1                                 (2.20)
%     th.kappa_het(...)         mean-reversion rate                                  (4.5)
%     th.stable_sign(theta,n,sigma)  +1 if the bump at x=0 is stable, else -1
%     th.Veff(...)              potential amplitude, printed (B.9) and derived
%     th.Deff(D,V)              effective diffusion  D / I0(2V/D)^2                  (4.6)
    th.Ahom   = @(theta,s) sqrt(1+theta) + s*sqrt(1-theta);
    th.ahom   = @(theta,s) acos(theta ./ (sqrt(1+theta) + s*sqrt(1-theta)));
    th.Dhom   = @(ep,theta) ep*pi ./ (2 + 2*sqrt(1-theta.^2));
    th.a_input     = @a_input;
    th.kappa_input = @(I0,n,a) n*I0.*sin(n*a) ./ (2*sin(a).^2 + n*I0.*sin(n*a));
    th.ou_var      = @(D,kap,t) D/(2*kap) * (1 - exp(-2*kap*t));
    th.Ahet        = @Ahet;
    th.a_het       = @a_het;
    th.lambda_odd  = @lambda_odd;
    th.lambda_odd_closed = @lambda_odd_closed;
    th.kappa_het   = @(theta,sg,n,sig) -th.lambda_odd(theta,sg,n,sig);
    th.stable_sign = @stable_sign;
    th.Veff        = @Veff;
    th.Deff        = @(D,V) D ./ besseli(0, 2*V./D).^2;
end

% ---------------------------------------------------------------------
function r = a_input(I0,n,theta)
% roots of sin(2a) + I0 cos(n a) = theta on (0,pi/2); r = [narrow wide].
% Stable (wide) branch has g'(a) < 0.
    g  = @(a) sin(2*a) + I0*cos(n*a) - theta;
    r  = roots_scan(g, 0, pi);
    r  = sort(r);
    if numel(r) == 1, r = [NaN r]; end
    if isempty(r),    r = [NaN NaN]; end
    r = [r(1) r(end)];
end

function A = Ahet(a,sg,n,sigma)
% A_+- = int_{-a}^{a} cos y (1 +- sigma cos(n y)) dy
    if n == 1
        Q = a + sin(a).*cos(a);
    else
        Q = sin((n+1)*a)/(n+1) + sin((n-1)*a)/(n-1);
    end
    A = 2*sin(a) + sg*sigma*Q;
end

function r = a_het(theta,sg,n,sigma)
% roots of A(a) cos(a) = theta, i.e. U(+-a) = theta; r = [narrow wide]
    g = @(a) Ahet(a,sg,n,sigma).*cos(a) - theta;
    r = sort(roots_scan(g, 0, pi/2));
    if numel(r) == 1, r = [NaN r]; end
    if isempty(r),    r = [NaN NaN]; end
    r = [r(1) r(end)];
end

function lam = lambda_odd(theta,sg,n,sigma)
% lambda_o = -1 + int sin^2 x (1 + sg sigma cos nx) f'(U) dx, Heaviside f'
% (2.18) with f'(U) = sum delta(x -+ a)/|U'(a)|,  U'(a) = -A sin a.
    r = a_het(theta,sg,n,sigma);  a = r(2);
    A = Ahet(a,sg,n,sigma);
    lam = -1 + 2*sin(a)*(1 + sg*sigma*cos(n*a)) / A;
end

function lam = lambda_odd_closed(theta,sg,n,sigma)
% closed form (2.20), n > 1
    r = a_het(theta,sg,n,sigma);  a = r(2);
    lam = sg*sigma*n*(n*sin(a)*cos(n*a) - cos(a)*sin(n*a)) / ...
        ((n^2-1)*sin(a) + sg*sigma*(n*cos(a)*sin(n*a) - sin(a)*cos(n*a)));
end

function sg = stable_sign(theta,n,sigma)
% Which bump class is the stable one: A_+ (centred at 0) or A_- (at pi/n)?
    if lambda_odd(theta,+1,n,sigma) < 0, sg = +1; else, sg = -1; end
end

function [Vprinted, Vderived] = Veff(theta,n,sigma)
% Amplitude of the cosine potential well V(D) = -V cos(n D).
%   Vprinted : formula (B.9) as printed in the paper.
%   Vderived : kappa/n^2 with kappa = -lambda_o (4.5), i.e. K(n D) = (kappa/n) sin(n D).
    sg = stable_sign(theta,n,sigma);
    r = a_het(theta,sg,n,sigma);  a = r(2);
    num = 2*sigma*(n*sin(a)*cos(n*a) - cos(a)*sin(n*a));
    den = n*(n^2-1)*sin(a) + sg*sigma*n*(n*cos(a)*sin(n*a) - sin(a)*cos(n*a));
    Vprinted = num/den;
    Vderived = -lambda_odd(theta,sg,n,sigma)/n^2;
end

% ---------------------------------------------------------------------
function r = roots_scan(g,lo,hi)
    a = linspace(lo+1e-9, hi-1e-9, 6001);
    v = g(a);
    idx = find(v(1:end-1).*v(2:end) < 0);
    r = zeros(1,numel(idx));
    for k = 1:numel(idx)
        r(k) = fzero(g, [a(idx(k)) a(idx(k)+1)]);
    end
    % tangential root (saddle-node) is dropped; report none/one honestly
end
