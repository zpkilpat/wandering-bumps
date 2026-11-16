function wb_save(fig,o,name)
%WB_SAVE  Save a figure as PNG and PDF in o.outdir.
    if ~o.save, return; end
    exportgraphics(fig, fullfile(o.outdir,[name '.png']), 'Resolution', 200);
    exportgraphics(fig, fullfile(o.outdir,[name '.pdf']), 'ContentType','vector');
end
