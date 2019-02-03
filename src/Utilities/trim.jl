# bound by a lower and upper quantiles
function trim!(x, percentage = 0.01)
    for j = 1:size(x,2)
        q = quantile(x[:,j], percentage)
        x[:,j] = max.(q,x[:,j])
        q = quantile(x[:,j], 1.0 - percentage)
        x[:,j] = min.(q,x[:,j])
    end
end    


