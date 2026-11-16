function figure7(varargin)
%FIGURE7  Fig. 7: bumps pinned by a stationary input I(x) = I0 cos(n x).
%   (a),(c) single realizations for n = 1, 2;  (b),(d) variance over M realizations
%   (red dashed) vs the Ornstein-Uhlenbeck theory (3.22) with kappa from (4.4)
%   (blue).  theta = 0.5, eps = 0.01, I0 = 0.1 (value inferred from the saturation
%   levels printed in the figure; the text does not state it).
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'src')));
    o = wb_opts(varargin{:});
    th = wb_theory();  theta = 0.5;  ep = 0.01;  I0 = 0.1;
    D = th.Dhom(ep,theta);
    fig = figure('Position',[100 100 900 760],'Color','w');
    Tv = wb_pick(o,[50 60],[30 30]);
    for k = 1:2
        n = k;
        out = nf_sim('eps',ep,'I0',I0,'nI',n,'T',100,'M',1,'keepfield',true,'seed',o.seed+k);
        ax = subplot(2,2,2*k-1); wb_field_panel(ax,out,100,sprintf('(%c)  n = %d',char('a'+2*k-2),n));
        s = nf_sim('eps',ep,'I0',I0,'nI',n,'T',Tv(k),'M',o.M,'seed',o.seed+10+k);
        v = mean(s.X.^2,1) - mean(s.X,1).^2;
        r = th.a_input(I0,n,theta);  a = r(2);
        kap = th.kappa_input(I0,n,a);
        ax = subplot(2,2,2*k); hold(ax,'on'); box(ax,'on');
        plot(ax,s.t,th.ou_var(D,kap,s.t),'b-','LineWidth',1.5);
        plot(ax,s.t,v,'r--','LineWidth',1.5);
        xlabel(ax,'t'); ylabel(ax,'\langle\Delta(t)^2\rangle'); set(ax,'FontSize',11);
        wb_label(ax,sprintf('(%c)',char('b'+2*k-2)));
        fprintf('figure7 n=%d: kappa = %.4f, saturation D/2kappa = %.4f\n',n,kap,D/(2*kap));
    end
    wb_save(fig,o,'fig7');
end
