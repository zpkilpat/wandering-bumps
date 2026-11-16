function figure5(varargin)
%FIGURE5  Fig. 5: wandering of bumps in the homogeneous ring.
%   (a),(b) single realizations, eps = 0.001 and 0.01 (theta = 0.5);
%   (c) <Delta^2> over M realizations vs the theory D(eps) t, eps = 0.01;
%   (d) D(eps) vs theta for eps = 0.001, 0.005, 0.01, theory (4.3) vs simulation.
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'src')));
    o = wb_opts(varargin{:});
    th = wb_theory();
    fig = figure('Position',[100 100 900 760],'Color','w');
    for k = 1:2
        ep = [0.001 0.01];
        out = nf_sim('eps',ep(k),'T',200,'M',1,'keepfield',true,'seed',o.seed+k);
        ax = subplot(2,2,k);  wb_field_panel(ax,out,200,sprintf('(%c)  \\epsilon = %g',char('a'+k-1),ep(k)));
    end
    % (c)
    T = wb_pick(o,50,30);
    out = nf_sim('eps',0.01,'T',T,'M',o.M,'seed',o.seed+10);
    v = mean(out.X.^2,1) - mean(out.X,1).^2;
    ax = subplot(2,2,3); hold(ax,'on'); box(ax,'on');
    plot(ax,out.t,th.Dhom(0.01,0.5)*out.t,'b-','LineWidth',1.5);
    plot(ax,out.t,v,'r--','LineWidth',1.5);
    xlabel(ax,'t'); ylabel(ax,'\langle\Delta(t)^2\rangle'); set(ax,'FontSize',11);
    wb_label(ax,'(c)'); text(ax,0.6*T,1.12*th.Dhom(0.01,0.5)*0.6*T,'D(\epsilon) t');
    % (d)
    thetas = wb_pick(o,0.05:0.1:0.95,[0.1 0.5 0.9]);
    eps_l  = [0.001 0.005 0.01];  Td = wb_pick(o,50,20);
    Dnum = NaN(numel(eps_l),numel(thetas));
    for i = 1:numel(eps_l)
        for j = 1:numel(thetas)
            s = nf_sim('eps',eps_l(i),'theta',thetas(j),'T',Td,'M',o.M,'seed',o.seed+100*i+j);
            vv = mean(s.X.^2,1) - mean(s.X,1).^2;
            Dnum(i,j) = sum(vv.*s.t)/sum(s.t.^2);              % slope through the origin
        end
    end
    ax = subplot(2,2,4); hold(ax,'on'); box(ax,'on');
    tt = linspace(0,0.995,300);
    for i = 1:numel(eps_l)
        plot(ax,tt,th.Dhom(eps_l(i),tt),'b-','LineWidth',1.4);
        plot(ax,thetas,Dnum(i,:),'ro','MarkerSize',5);
        text(ax,0.55,th.Dhom(eps_l(i),0.5)*1.25,sprintf('\\epsilon = %g',eps_l(i)));
    end
    xlabel(ax,'\theta'); ylabel(ax,'D(\epsilon)'); ylim(ax,[0 0.014]); xlim(ax,[0 1]);
    wb_label(ax,'(d)'); set(ax,'FontSize',11);
    wb_save(fig,o,'fig5');
end
