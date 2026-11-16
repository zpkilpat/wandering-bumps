function v = wb_pick(o,fullval,quickval)
%WB_PICK  Choose the paper-scale value or the quick-mode value.
    if o.quick, v = quickval; else, v = fullval; end
end
