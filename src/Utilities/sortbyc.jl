# sort by column
function sortbyc(x,col)
    x = sortslices(x,dims=1,by=x->x[col])
end

function sortbyc()
    a = [1 2 3; 3 2 1]
    sortbyc(a,3)[1,1] == 3
end    

