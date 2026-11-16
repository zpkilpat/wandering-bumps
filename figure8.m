function figure8(varargin)
%FIGURE8  Fig. 8: escape of a pinned bump between the two attractors of I0 cos(2x).
%   (a) one realization (a run with a hop to x ~ pi is selected by scanning seeds);
%   (b) mean switching time vs I0, fit to b exp(gamma I0) (paper: b = 750, gamma = 30).
%   A switch is recorded when the bump centre first comes within 0.1 pi of the other
%   attractor (|centre| >= 0.9 pi); the paper does not state its criterion.
%   These runs are long, so dt = 0.05 is used here (the paper's dt = 0.01 elsewhere).
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'src')));
    o = wb_opts(varargin{:});
    th = wb_theory(); %#ok<NASGU>
    dt = 0.05;  thr = 0.9*pi;  ep = 0.01;
    fig = figure('Position',[100 100 900 340],'Color','w');
    I0a  = wb_pick(o,0.05,0.03);  Ta = wb_pick(o,10000,4000);
    win  = wb_pick(o,[3000 7000],[800 3000]);
    chosen = [];
    for sd = 1:60
        out = nf_sim('eps',ep,'I0',I0a,'nI',2,'T',Ta,'dt',dt,'rec',40,'M',1,'keepfield',true,'seed',sd);
        i = find(abs(out.X(1,:)) >= thr, 1);
        if ~isempty(i) && out.t(i) > win(1) && out.t(i) < win(2), chosen = out; break; end
    end
    if isempty(chosen), chosen = out; warning('figure8: no hop in the window'); end
    ax = subplot(1,2,1); wb_field_panel(ax,chosen,Ta,'(a)');

    I0s = wb_pick(o,linspace(0,0.1,13),linspace(0,0.05,5));
    Mw  = wb_pick(o,20,6);  Tcap = wb_pick(o,60000,8000);
    Ts  = NaN(size(I0s));
    for j = 1:numel(I0s)
        s = nf_sim('eps',ep,'I0',I0s(j),'nI',2,'T',Tcap,'dt',dt,'rec',40,'M',Mw,'stop','switch', ...
                   'stopthr',thr,'seed',o.seed+j);
        e = s.event;  e(isnan(e)) = Tcap;  Ts(j) = mean(e);
        fprintf('figure8: I0 = %.4f  mean T_switch = %.0f  (%d/%d switched)\n',I0s(j),Ts(j),sum(~isnan(s.event)),Mw);
    end
    p = fminsearch(@(p) sum((p(1)*exp(p(2)*I0s) - Ts).^2),[750 30]);
    ax = subplot(1,2,2); hold(ax,'on'); box(ax,'on');
    tt = linspace(0,max(I0s),200);
    plot(ax,tt,p(1)*exp(p(2)*tt),'b-','LineWidth',1.5);
    plot(ax,I0s,Ts,'ro','MarkerSize',5);
    xlabel(ax,'I_0'); ylabel(ax,'T_{switch}'); set(ax,'FontSize',11);
    wb_note(ax,sprintf('fit: b = %.0f, \\gamma = %.1f',p(1),p(2)),'br');
    wb_label(ax,'(b)');
    fprintf('figure8: fit b = %.0f, gamma = %.1f  (paper: 750, 30)\n',p(1),p(2));
    wb_save(fig,o,'fig8');
end
