function figure1(varargin)
%FIGURE1  Fig. 1: input-locked bumps in the deterministic ring model.
%   (a) bump half-width a vs input strength I0 for n=1 (2.7) and n=2 (2.8)-(2.9);
%   stable (wide, black) and unstable (narrow, grey) branches.  theta = 0.5.
%   (b) stable / unstable bumps for the bimodal input, I0 = 0.4.
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'src')));
    o = wb_opts(varargin{:});
    th = wb_theory();  theta = 0.5;
    I0 = 0:0.005:0.5;
    aN = NaN(numel(I0),2);  aW = aN;
    for n = 1:2
        for k = 1:numel(I0)
            r = th.a_input(I0(k), n, theta);  aN(k,n) = r(1);  aW(k,n) = r(2);
        end
    end
    fig = figure('Position',[100 100 900 340],'Color','w');
    ax = subplot(1,2,1); hold(ax,'on'); box(ax,'on');
    plot(ax,I0,aW(:,1),'k-','LineWidth',1.5);  plot(ax,I0,aW(:,2),'k-','LineWidth',1.5);
    plot(ax,I0,aN(:,1),'-','Color',[.5 .5 .5],'LineWidth',1.5);
    plot(ax,I0,aN(:,2),'-','Color',[.5 .5 .5],'LineWidth',1.5);
    text(ax,0.33,1.33,'n = 1'); text(ax,0.33,1.04,'n = 2'); text(ax,0.3,0.2,'n = 1,2');
    xlabel(ax,'I_0'); ylabel(ax,'a'); ylim(ax,[0 1.45]); wb_label(ax,'(a)');
    set(ax,'FontSize',11);

    ax = subplot(1,2,2); hold(ax,'on'); box(ax,'on');
    I = 0.4;  r = th.a_input(I,2,theta);  x = linspace(-pi,pi,800);
    for s = [1 -1]
        plot(ax,x, s*2*sin(r(2))*cos(x) + I*cos(2*x),'k-','LineWidth',1.5);
        plot(ax,x, s*2*sin(r(1))*cos(x) + I*cos(2*x),'-','Color',[.6 .6 .6],'LineWidth',1.5);
    end
    plot(ax,x,I*cos(2*x),'k:','LineWidth',1.2);
    plot(ax,[-pi pi],[theta theta],'-','Color',[.8 .8 .8]);
    text(ax,0.3,2.05,'stable'); text(ax,0.3,0.95,'unstable');
    xlabel(ax,'x'); ylabel(ax,'U(x)'); xlim(ax,[-pi pi]); ylim(ax,[-1.6 2.2]);
    wb_label(ax,'(b)'); set(ax,'FontSize',11);
    wb_save(fig,o,'fig1');
end
