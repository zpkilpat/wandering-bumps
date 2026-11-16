function figure10(varargin)
%FIGURE10  Fig. 10: reduced effective diffusion, high-frequency modulation.
%   w = (1 + sigma cos(n y)) cos(x - y); theta = 0.5, sigma = 0.1, eps = 0.01.
%   (a) one realization, n = 8, attractors in cyan;   (b) variance vs D_eff t;
%   (c) density P(Delta,t) at t = 400 over 5000 realizations;
%   (d) D_eff vs n: theory (4.6) (blue), simulation (red dash-dot), D(eps) (black).
%
%   Theory: D_eff = D / I_0(2V/D)^2 with a cosine well V(Delta) = -V cos(n Delta).
%   V = kappa/n^2, kappa = -lambda_o (4.5), which reproduces the paper's blue curve.
%   The printed (B.9) is -2 kappa/n^2 (wrong sign and a factor 2); see README.
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'src')));
    o = wb_opts(varargin{:});
    th = wb_theory();  theta = 0.5;  ep = 0.01;  sigma = 0.1;
    D = th.Dhom(ep,theta);
    n = 8;  sg = th.stable_sign(theta,n,sigma);
    [~,V] = th.Veff(theta,n,sigma);  Deff8 = th.Deff(D,V);
    fig = figure('Position',[100 100 900 760],'Color','w');
    cy = [0 .8 .9];
    if sg > 0, cen = (2*pi/n)*(-n/2+1:n/2); else, cen = (2*pi/n)*((-n/2:n/2-1)+0.5); end

    out = nf_sim('eps',ep,'n',n,'sigma',sigma,'T',400,'M',1,'keepfield',true,'seed',o.seed);
    ax = subplot(2,2,1); wb_field_panel(ax,out,400,'(a)  n = 8'); hold(ax,'on');
    for c = cen, plot(ax,[c c],[0 400],'-','Color',cy,'LineWidth',1); end

    Mbig = wb_pick(o,5000,500);  Tbig = wb_pick(o,400,100);
    s = nf_sim('eps',ep,'n',n,'sigma',sigma,'T',Tbig,'M',Mbig,'seed',o.seed+1,'rec',50);
    v = mean(s.X.^2,1) - mean(s.X,1).^2;
    ax = subplot(2,2,2); hold(ax,'on'); box(ax,'on');
    m = s.t <= 30;
    plot(ax,s.t(m),Deff8*s.t(m),'b-','LineWidth',1.5);
    plot(ax,s.t(m),v(m),'r--','LineWidth',1.5);
    xlabel(ax,'t'); ylabel(ax,'\langle\Delta^2\rangle'); set(ax,'FontSize',11);
    wb_label(ax,'(b)');

    ax = subplot(2,2,3); hold(ax,'on'); box(ax,'on');
    edges = linspace(-pi,pi,61);
    pos = angle(exp(1i*(s.X(:,end))));
    h = histcounts(pos,edges)/numel(pos);
    bar(ax,(edges(1:end-1)+edges(2:end))/2,h,1,'FaceColor','k','EdgeColor','k');
    for c = cen, plot(ax,[c c],[0 max(h)*1.1],'-','Color',cy,'LineWidth',1); end
    xlim(ax,[-pi pi]); xlabel(ax,'\Delta'); ylabel(ax,'P(\Delta,t)'); set(ax,'FontSize',11);
    wb_label(ax,sprintf('(c)  t = %g',Tbig));

    ns = 2:10;  Td = wb_pick(o,200,60);  Md = wb_pick(o,1000,100);
    Dth = NaN(size(ns));  Dnum = NaN(size(ns));
    for j = 1:numel(ns)
        nn = ns(j);
        [~,Vn] = th.Veff(theta,nn,sigma);  Dth(j) = th.Deff(D,Vn);
        r = nf_sim('eps',ep,'n',nn,'sigma',sigma,'T',Td,'M',Md,'seed',o.seed+20+j,'rec',50);
        vv = mean(r.X.^2,1) - mean(r.X,1).^2;
        w = r.t >= Td/2;
        pf = polyfit(r.t(w),vv(w),1);  Dnum(j) = pf(1);        % long-time slope = D_eff
    end
    ax = subplot(2,2,4); hold(ax,'on'); box(ax,'on');
    h1 = plot(ax,[2 10],[D D],'k-');
    h2 = plot(ax,ns,Dth,'b-','LineWidth',1.5);
    h3 = plot(ax,ns,Dnum,'r.-.','MarkerSize',14);
    legend(ax,[h1 h2 h3],{'D(\epsilon)','theory (4.6)','simulation'},'Location','southeast','Box','off');
    xlabel(ax,'n'); ylabel(ax,'D_{eff}'); xlim(ax,[2 10]); ylim(ax,[0 0.01]); set(ax,'FontSize',11);
    wb_label(ax,'(d)');
    fprintf('figure10: n  Dth  Dnum\n'); fprintf('   %2d  %.5f  %.5f\n',[ns; Dth; Dnum]);
    wb_save(fig,o,'fig10');
end
