function figure2(varargin)
%FIGURE2  Fig. 2: 2n bump locations with synaptic modulation n.
%   (a) bump centres x = m pi/n; stable (blue filled) and unstable (red circles).
%   (b) stable (A_+, blue) and unstable (A_-, red dashed) bumps for n = 3.
%   theta = 0.5, sigma = 0.2.
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'src')));
    o = wb_opts(varargin{:});
    th = wb_theory();  theta = 0.5;  sigma = 0.2;
    fig = figure('Position',[100 100 900 340],'Color','w');
    ax = subplot(1,2,1); hold(ax,'on'); box(ax,'on');
    for n = 1:6
        plot(ax,[-pi pi],[n n],'k-');
        for m = -n+1:n
            sg = 1 - 2*mod(m,2);                         % even m: A_+, odd m: A_-
            stable = th.lambda_odd(theta,sg,n,sigma) < 0;
            if stable
                plot(ax,m*pi/n,n,'o','MarkerFaceColor',[.1 .2 .7],'MarkerEdgeColor',[.1 .2 .7],'MarkerSize',6);
            else
                plot(ax,m*pi/n,n,'o','MarkerFaceColor','w','MarkerEdgeColor',[.8 .1 .1],'MarkerSize',6);
            end
        end
    end
    plot(ax,[-pi pi],[0 0],'k-');
    xlim(ax,[-3.4 3.4]); ylim(ax,[-0.3 6.5]); yticks(ax,0:6); xticks(ax,-3:3);
    xlabel(ax,'x'); ylabel(ax,'n'); wb_label(ax,'(a)'); set(ax,'FontSize',11);

    ax = subplot(1,2,2); hold(ax,'on'); box(ax,'on');
    n = 3;  x = linspace(-pi,pi,1000);
    aP = th.a_het(theta,+1,n,sigma);  AP = th.Ahet(aP(2),+1,n,sigma);
    aM = th.a_het(theta,-1,n,sigma);  AM = th.Ahet(aM(2),-1,n,sigma);
    for m = -n+1:n
        c = m*pi/n;
        if mod(m,2) == 0
            plot(ax,x,AP*cos(x-c),'-','Color',[.1 .2 .7],'LineWidth',1.3);
        else
            plot(ax,x,AM*cos(x-c),'--','Color',[.8 .1 .1],'LineWidth',1.3);
        end
    end
    xlim(ax,[-pi pi]); ylim(ax,[-2.2 2.2]); xticks(ax,-3:3);
    xlabel(ax,'x'); ylabel(ax,'U(x)'); wb_label(ax,'(b)'); set(ax,'FontSize',11);
    wb_save(fig,o,'fig2');
end
