using Statistics
function BoundByQuantiles!(data, margin=0.005)
for j = 1:size(data,2)
    q = quantile(data[:,j],margin)
    data[:,j] = max.(q,data[:,j])
    q = quantile(data[:,j],1.0-margin)
    data[:,j] = min.(q,data[:,j])
end
end

