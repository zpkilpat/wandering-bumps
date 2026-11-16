function wb_label(ax,str,col)
%WB_LABEL  Panel label ('(a)', ...) in the top-left corner of an axes.
%   Call it after the panel is drawn and its limits are set: the position is
%   computed in data units from the current limits, which are then frozen.
    if nargin < 3, col = 'k'; end
    xl = xlim(ax);  yl = ylim(ax);  xlim(ax,xl);  ylim(ax,yl);
    text(ax, xl(1)+0.03*diff(xl), yl(2)-0.03*diff(yl), str, 'Color',col, ...
         'FontWeight','bold','FontSize',11,'VerticalAlignment','top', ...
         'HorizontalAlignment','left','Clipping','off');
end
