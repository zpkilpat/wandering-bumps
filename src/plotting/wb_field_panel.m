function wb_field_panel(ax,out,tmax,label)
%WB_FIELD_PANEL  Space-time plot of U(x,t) with the bump centre superimposed
%   (the layout used in Figs 5, 6, 7, 8, 9, 10 of the paper).
    imagesc(ax, out.x, out.t, out.field'); set(ax,'YDir','normal');
    colormap(ax, hot); clim(ax,[0 2]); hold(ax,'on');
    c = angle(exp(1i*(out.c0 + out.X(1,:))));          % centre wrapped to (-pi,pi]
    c([false abs(diff(c)) > pi]) = NaN;                 % break the line at the seam
    if ~isnan(out.event(1)) && ~strcmp(out.params.stop,'switch')
        c(out.t >= out.event(1)) = NaN;                 % bump is extinct: stop the line
    end
    plot(ax, c, out.t, '-', 'Color',[0.55 0.1 0.8], 'LineWidth',1.5);
    xlim(ax,[-pi pi]); ylim(ax,[0 tmax]); xticks(ax,-3:3);
    xlabel(ax,'x'); ylabel(ax,'t'); set(ax,'FontSize',11);
    wb_label(ax,label,'w');
end
