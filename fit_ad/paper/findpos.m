function [idx,weight]=findpos(x_axis,x_current)

%%

if x_current<x_axis(1)
    idx=[1 2];
    weight=[ 0 0];
    return
end


if x_current>=x_axis(end)
    idx=[length(x_axis)-1 length(x_axis)];
    weight=[ 0 0];
    return
end

idx_end=find(x_current<x_axis,1);


idx=[idx_end-1 idx_end];

L_el=x_axis(idx(2))-x_axis(idx(1));


weight=[ 1-(x_current-x_axis(idx_end-1))/L_el (x_current-x_axis(idx_end-1))/L_el];

