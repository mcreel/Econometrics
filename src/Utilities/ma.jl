# compute moving average using p most recent values, including current value
function ma(x, p)
    m = similar(x)
    for i = p:size(x,1)
        m[i] = mean(x[i-p+1:i])
    end
    return m
end


