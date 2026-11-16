function figure3(varargin)
%FIGURE3  Fig. 3: bumps with heterogeneous weights, modulation n = 1.
%   (a) half-width a vs heterogeneity amplitude sigma, from (2.16);
%   (b) profiles for sigma = 0.2.  theta = 0.5.
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'src')));
    o = wb_opts(varargin{:});
    th = wb_theory();  theta = 0.5;
    sg_ = 0:0.002:0.6;
    aP = NaN(numel(sg_),2);  aM = aP;
    for k = 1:numel(sg_)
        aP(k,:) = th.a_het(theta,+1,1,sg_(k));
        aM(k,:) = th.a_het(theta,-1,1,sg_(k));
    end
    fig = figure('Position',[100 100 900 340],'Color','w');
    ax = subplot(1,2,1); hold(ax,'on'); box(ax,'on');
    plot(ax,sg_,aP(:,2),'-','Color',[.1 .2 .7],'LineWidth',1.6);       % x=0 wide, stable
    plot(ax,sg_,aP(:,1),'-','Color',[.5 .5 .5],'LineWidth',1.6);       % x=0 narrow
    plot(ax,sg_,aM(:,2),'--','Color',[.8 .1 .1],'LineWidth',1.6);      % x=pi wide
    plot(ax,sg_,aM(:,1),':','Color',[.5 .5 .5],'LineWidth',1.6);       % x=pi narrow
    text(ax,0.3,1.42,'x = 0'); text(ax,0.38,0.92,'x = \pi'); text(ax,0.3,0.12,'x = 0');
    xlabel(ax,'\sigma'); ylabel(ax,'a'); xlim(ax,[0 0.6]); ylim(ax,[0 1.65]);
    wb_label(ax,'(a)'); set(ax,'FontSize',11);

    ax = subplot(1,2,2); hold(ax,'on'); box(ax,'on');
    sigma = 0.2;  x = linspace(-pi,pi,1000);
    rP = th.a_het(theta,+1,1,sigma);  rM = th.a_het(theta,-1,1,sigma);
    plot(ax,x, th.Ahet(rP(2),+1,1,sigma)*cos(x),'-','Color',[.1 .2 .7],'LineWidth',1.6);
    plot(ax,x, th.Ahet(rP(1),+1,1,sigma)*cos(x),'-','Color',[.5 .5 .5],'LineWidth',1.6);
    plot(ax,x,-th.Ahet(rM(2),-1,1,sigma)*cos(x),'--','Color',[.8 .1 .1],'LineWidth',1.6);
    plot(ax,x,-th.Ahet(rM(1),-1,1,sigma)*cos(x),'--','Color',[.5 .5 .5],'LineWidth',1.6);
    xlim(ax,[-pi pi]); ylim(ax,[-2.2 2.4]); xticks(ax,-3:3);
    xlabel(ax,'x'); ylabel(ax,'U(x)'); wb_label(ax,'(b)'); set(ax,'FontSize',11);
    wb_save(fig,o,'fig3');
end
