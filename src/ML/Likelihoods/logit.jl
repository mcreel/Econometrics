function logit(theta, y, x)
    p = 1.0./(1.0 .+ exp.(-x*theta))
    obj = y.*log.(p) .+ (log.(1.0 .- p)).*(1.0 .- y)
end    
