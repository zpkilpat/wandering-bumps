function figure9(varargin)
%FIGURE9  Fig. 9: pinning by low-frequency synaptic heterogeneity.
%   w = (1 + sigma cos(n y)) cos(x - y), n = 2 and 3, sigma = 0.1, eps = 0.01, theta = 0.5.
%   (a),(c) single realizations;  (b),(d) variance over M realizations (red dashed)
%   vs the Ornstein-Uhlenbeck theory (3.22) with kappa from (4.5) (blue).
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'src')));
    o = wb_opts(varargin{:});
    th = wb_theory();  theta = 0.5;  ep = 0.01;  sigma = 0.1;
    D = th.Dhom(ep,theta);
    fig = figure('Position',[100 100 900 760],'Color','w');
    for k = 1:2
        n = k + 1;
        out = nf_sim('eps',ep,'n',n,'sigma',sigma,'T',100,'M',1,'keepfield',true,'seed',o.seed+k);
        ax = subplot(2,2,2*k-1); wb_field_panel(ax,out,100,sprintf('(%c)  n = %d',char('a'+2*k-2),n));
        s = nf_sim('eps',ep,'n',n,'sigma',sigma,'T',30,'M',o.M,'seed',o.seed+10+k);
        v = mean(s.X.^2,1) - mean(s.X,1).^2;
        kap = th.kappa_het(theta,+1,n,sigma);
        ax = subplot(2,2,2*k); hold(ax,'on'); box(ax,'on');
        plot(ax,s.t,th.ou_var(D,kap,s.t),'b-','LineWidth',1.5);
        plot(ax,s.t,v,'r--','LineWidth',1.5);
        xlabel(ax,'t'); ylabel(ax,'\langle\Delta(t)^2\rangle'); set(ax,'FontSize',11);
        wb_label(ax,sprintf('(%c)',char('b'+2*k-2)));
        fprintf('figure9 n=%d: kappa = %.4f, saturation D/2kappa = %.4f\n',n,kap,D/(2*kap));
    end
    wb_save(fig,o,'fig9');
end
