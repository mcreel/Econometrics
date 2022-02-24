#  Negative Binomial type 1 and type II log likelihood
#  The parameterization follows
#  Deb and Trivedi, J. Appl. Econometrics,
#  V. 12, 1997, pp. 313 - 336.
using SpecialFunctions, Distributions
function  negbin(θ, y, x, nbtype)
    n, k = size(x)
    β = θ[1:k]
    eps = 1e-8
    λ = exp.(eps .+ x*β)
    α = exp(eps .+ θ[end])
    nbtype == 1 ? r = λ/α : r = ones(n)/α
    p = r ./ (r + λ)
    logdensity = log.(pdf.(NegativeBinomial.(r, p),y))    
end


