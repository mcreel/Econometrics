function logit(θ, y, x)
    p = 1.0./(1.0 .+ exp.(-x*θ))
    y.*log.(p) .+ (log.(1.0 .- p)).*(1.0 .- y)
end    
