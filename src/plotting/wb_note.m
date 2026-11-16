function wb_note(ax,str,loc)
%WB_NOTE  Small note inside an axes: loc = 'tr' (top right) or 'br' (bottom right).
%   Same conventions as WB_LABEL (call after the limits are set).
    if nargin < 3, loc = 'tr'; end
    xl = xlim(ax);  yl = ylim(ax);  xlim(ax,xl);  ylim(ax,yl);
    x = xl(2) - 0.03*diff(xl);
    if strcmp(loc,'br'), y = yl(1) + 0.04*diff(yl); va = 'bottom';
    else,                y = yl(2) - 0.03*diff(yl); va = 'top'; end
    text(ax,x,y,str,'FontSize',10,'VerticalAlignment',va,'HorizontalAlignment','right');
end
