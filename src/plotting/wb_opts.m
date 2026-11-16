function o = wb_opts(varargin)
%WB_OPTS  Options shared by figure1 ... figure10 (name-value pairs).
%   'quick' (false)  small realization counts and short runs, a few minutes in total
%   'fast'  (false)  use the exact 2D reduction (nf_sim_fast) wherever there is no input
%   'M'     (1000)   realizations for variance estimates (paper: 1000)
%   'save'  (true)   write PNG and PDF to 'outdir'
%   'outdir'         default: <repo>/output/figs
%   'seed'  (1)      base random seed
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    o = struct('quick',false,'fast',false,'M',1000,'save',true, ...
               'outdir',fullfile(root,'output','figs'),'seed',1);
    for k = 1:2:numel(varargin), o.(varargin{k}) = varargin{k+1}; end
    if o.quick && ~any(strcmp(varargin(1:2:end),'M')), o.M = 100; end
    if ~exist(o.outdir,'dir'), mkdir(o.outdir); end
    setappdata(0,'wb_fast',logical(o.fast));          % read by nf_sim
end
