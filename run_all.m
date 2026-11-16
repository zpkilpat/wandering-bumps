function run_all(varargin)
%RUN_ALL  Reproduce every figure of
%     Z. P. Kilpatrick and B. Ermentrout, "Wandering bumps in stochastic neural fields",
%     SIAM J. Appl. Dyn. Syst. 12(1):61-94, 2013.
%
%   run_all                    paper-scale runs with the grid solver (about 2 hours)
%   run_all('fast',true)       exact 2D reduction where there is no input (about 10 minutes,
%                              mostly Figs. 7 and 8, which have an input)
%   run_all('quick',true)      fewer realizations and shorter runs (a few minutes)
%   run_all('figs',[5 7])      only some figures
%   Other options: 'outdir', 'M', 'seed' (see src/plotting/wb_opts.m).
    p = inputParser;  p.KeepUnmatched = true;
    addParameter(p,'figs',1:10);
    parse(p,varargin{:});
    opts = varargin;
    k = find(strcmp(opts,'figs'));  if ~isempty(k), opts(k:k+1) = []; end
    for k = p.Results.figs
        tic;  fprintf('--- Figure %d\n',k);
        feval(sprintf('figure%d',k),opts{:});
        fprintf('    done in %.0f s\n',toc);
    end
end
