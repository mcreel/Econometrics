function dstats(x, rnames="";short=false)
    k = size(x,2)
    if rnames==""
        rnames = 1:k
        rnames = rnames'
    end
    m = mean(x,dims=1)
    mm = median(x,dims=1)
    s = std(x,dims=1)
    sk = m-m
    kt = m-m
    mn = minimum(x,dims=1)
    mx = maximum(x,dims=1)
    q05 = fill(0.0,k)
    q95 = fill(0.0,k)
    if short == false
        for i = 1:size(x,2) q05[i] = quantile(x[:,i], 0.05) end
        for i = 1:size(x,2) q95[i] = quantile(x[:,i], 0.95) end
        cnames = ["  mean", " median","  std", "min", "max", "q05", "q95"]
        stats = [m' mm' s' mn' mx' q05 q95] 
        prettyprint(stats, cnames, rnames)
    else
        cnames = ["  mean", " median", "  std", "min", "max"]
        stats = [m' mm' s' mn' mx'] 
        prettyprint(stats, cnames, rnames)
    end
    return stats
end

# test dstats
function dstats()
    a = [1;2;3]
    b = dstats(a)[:];
    c = [2.0;2.0;1.0;1.0;3.0;1.1;2.9]
    b == c
end    



