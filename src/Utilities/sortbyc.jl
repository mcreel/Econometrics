# sort by column
function sortbyc(x,col)
    x = sortslices(x,dims=1,by=x->x[col])
end
