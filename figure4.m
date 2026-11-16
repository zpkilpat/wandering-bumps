function figure4(varargin)
%FIGURE4  Fig. 4: eigenvalue lambda_o of odd (shift) perturbations of the
%   bump at x = 0 in the network with synaptic modulation (1.4).  theta = 0.5.
%   (a) lambda_o vs sigma for n = 1, 2, 4.   (b) lambda_o vs n for sigma = 0.2.
%
%   NOTE on n = 1.  The paper's figure for n = 1 follows the printed closed form
%   (2.19),  -sigma a / (2 sin a + sigma a + sigma sin 2a).  Evaluating (2.18)
%   directly (and measuring the decay of a shifted bump in a deterministic
%   simulation) gives  sigma (sin a cos a - a) / A,  A = 2 sin a + sigma (a + sin a cos a),
%   which is slightly less negative.  Both are drawn: solid = as printed (matches
%   the paper's figure), dashed = direct evaluation of (2.18).
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'src')));
    o = wb_opts(varargin{:});
    th = wb_theory();  theta = 0.5;
    lam1 = @(sig) printed219(th,theta,sig);
    fig = figure('Position',[100 100 900 340],'Color','w');
    ax = subplot(1,2,1); hold(ax,'on'); box(ax,'on');
    sg_ = 0:0.01:0.6;
    L1 = arrayfun(lam1,sg_);
    L1d = arrayfun(@(s) th.lambda_odd(theta,1,1,s), sg_);
    L2 = arrayfun(@(s) th.lambda_odd(theta,1,2,s), sg_);
    L4 = arrayfun(@(s) th.lambda_odd(theta,1,4,s), sg_);
    plot(ax,sg_,L1,'-','Color',[.1 .2 .7],'LineWidth',1.5);
    plot(ax,sg_,L1d,'--','Color',[.1 .2 .7],'LineWidth',1);
    plot(ax,sg_,L2,'-','Color',[.1 .2 .7],'LineWidth',1.5);
    plot(ax,sg_,L4,'-','Color',[.8 .1 .1],'LineWidth',1.5);
    plot(ax,[0 0.6],[0 0],'-','Color',[.8 .8 .8]);
    text(ax,0.38,-0.17,'n = 1'); text(ax,0.25,-0.3,'n = 2'); text(ax,0.4,0.18,'n = 4');
    xlabel(ax,'\sigma'); ylabel(ax,'\lambda_o'); xlim(ax,[0 0.6]); ylim(ax,[-0.5 0.32]);
    wb_label(ax,'(a)'); set(ax,'FontSize',11);

    ax = subplot(1,2,2); hold(ax,'on'); box(ax,'on');
    sigma = 0.2;  ns = 1:12;
    lam = arrayfun(@(n) th.lambda_odd(theta,1,n,sigma), ns);
    lam(1) = lam1(sigma);                                        % as printed (2.19)
    plot(ax,ns,lam,'-','Color',[.1 .2 .7],'LineWidth',1.3);
    plot(ax,ns,lam,'o','MarkerFaceColor',[.6 .85 .1],'MarkerEdgeColor',[.3 .5 0],'MarkerSize',6);
    plot(ax,1,th.lambda_odd(theta,1,1,sigma),'x','Color',[.1 .2 .7]);   % direct (2.18)
    plot(ax,[0 13],[0 0],'-','Color',[.8 .8 .8]);
    text(ax,5,0.27,'x = 0 : unstable'); text(ax,5,-0.27,'x = 0 : stable');
    xlabel(ax,'n'); ylabel(ax,'\lambda_o'); xlim(ax,[0.5 12.5]); ylim(ax,[-0.3 0.32]);
    wb_label(ax,'(b)'); set(ax,'FontSize',11);
    wb_save(fig,o,'fig4');
end

function lam = printed219(th,theta,sigma)
    r = th.a_het(theta,+1,1,sigma);  a = r(2);
    lam = -sigma*a / (2*sin(a) + sigma*a + sigma*sin(2*a));      % (2.19) as printed
end
