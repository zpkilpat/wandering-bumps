function figure6(varargin)
%FIGURE6  Fig. 6: extinction of bumps near the saddle node (theta_SN = 1).
%   (a) one realization, theta = 0.95, eps = 0.01 (a run that goes extinct near t ~ 65
%       is selected by scanning seeds);
%   (b) mean extinction time vs theta over M realizations, least-squares fit to
%       b exp(gamma |theta - theta_SN|) (paper: b ~ 10, gamma ~ 33).
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'src')));
    o = wb_opts(varargin{:});
    fig = figure('Position',[100 100 900 340],'Color','w');
    chosen = [];
    for sd = 1:300
        out = nf_sim('theta',0.95,'eps',0.01,'T',70,'M',1,'keepfield',true,'seed',sd);
        if ~isnan(out.event) && out.event > 55 && out.event < 75, chosen = out; break; end
    end
    if isempty(chosen), chosen = out; warning('figure6: no run went extinct near t=65'); end
    ax = subplot(1,2,1); wb_field_panel(ax,chosen,70,'(a)');

    thetas = wb_pick(o,0.925:0.005:0.995,0.94:0.02:0.99);
    Tcap   = wb_pick(o,3000,600);
    Tm = NaN(size(thetas));  ncens = zeros(size(thetas));
    for j = 1:numel(thetas)
        s = nf_sim('theta',thetas(j),'eps',0.01,'T',Tcap,'M',o.M,'stop','extinct','seed',o.seed+j,'rec',20);
        e = s.event;  ncens(j) = sum(isnan(e));  e(isnan(e)) = Tcap;
        Tm(j) = mean(e);
    end
    f = @(p) sum((p(1)*exp(p(2)*abs(thetas-1)) - Tm).^2);
    p = fminsearch(f,[10 33]);
    ax = subplot(1,2,2); hold(ax,'on'); box(ax,'on');
    tt = linspace(min(thetas),1,200);
    plot(ax,tt,p(1)*exp(p(2)*abs(tt-1)),'b-','LineWidth',1.5);
    plot(ax,thetas,Tm,'ro','MarkerSize',5);
    plot(ax,[1 1],[0 max(Tm)*1.1],'k-');
    text(ax,0.99,max(Tm)*0.5,'saddle-node','Rotation',90,'HorizontalAlignment','center');
    xlabel(ax,'\theta'); ylabel(ax,'T_{extinct}'); xlim(ax,[min(thetas)-0.002 1.003]);
    wb_label(ax,'(b)'); set(ax,'FontSize',11);
    wb_note(ax,sprintf('fit: b = %.1f, \\gamma = %.1f',p(1),p(2)),'tr');
    fprintf('figure6: fit b = %.2f, gamma = %.2f  (paper: ~10, ~33); censored runs: %d\n',p(1),p(2),sum(ncens));
    wb_save(fig,o,'fig6');
end
