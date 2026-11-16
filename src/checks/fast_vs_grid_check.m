function fast_vs_grid_check(M)
%FAST_VS_GRID_CHECK  Compare the 2D reduced solver (nf_sim_fast) with the grid solver (nf_sim).
%   Prints the diffusion coefficient (free bump), the OU saturation (n = 2, 3 pinning)
%   and run times; both should agree with each other and with theory to within
%   sampling error ~ sqrt(2/M).
    if nargin < 1, M = 1000; end
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'..')));
    th = wb_theory();  ep = 0.01;  D = th.Dhom(ep,0.5);
    fprintf('theory D = %.5f\n',D);
    cases = {'free (sigma=0), slope/D', 0, 1, 50; ...
             'n=2 sigma=.1, plateau/theory', 0.1, 2, 30; ...
             'n=3 sigma=.1, plateau/theory', 0.1, 3, 30};
    for c = 1:size(cases,1)
        [lab,sg,n,T] = cases{c,:};
        for fast = [false true]
            tic;
            if fast, o = nf_sim_fast('eps',ep,'sigma',sg,'n',n,'T',T,'M',M,'seed',7);
            else,    o = nf_sim('eps',ep,'sigma',sg,'n',n,'T',T,'M',M,'seed',7); end
            el = toc;
            v = mean(o.X.^2,1) - mean(o.X,1).^2;
            if sg == 0
                val = (sum(v.*o.t)/sum(o.t.^2))/D;
            else
                kap = th.kappa_het(0.5,+1,n,sg);
                val = mean(v(o.t > T/2)) / (D/(2*kap)*mean(1-exp(-2*kap*o.t(o.t > T/2))));
            end
            fprintf('%-32s %-5s  %.3f   (%.1f s)\n',lab,ternary(fast,'fast','grid'),val,el);
        end
    end
end
function s = ternary(c,a,b), if c, s = a; else, s = b; end, end
