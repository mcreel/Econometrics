#  Negative Binomial type 1 and type II log likelihood
#  The parameterization follows
#  Deb and Trivedi, J. Appl. Econometrics,
#  V. 12, 1997, pp. 313 - 336.
using SpecialFunctions, Distributions
function  negbin(θ, y, x, nbtype)
    n, k = size(x)
    β = θ[1:k]
    λ = exp.(x*β)
    α = exp(θ[end])
    nbtype == 1 ? r = α./λ : r = ones(n)/α
    p = r ./ (r + λ)
    ((r > zero(r)) & (p > zero(p)) & ((1.0 .-p) > zero(p))) ? log.(pdf.(NegativeBinomial.(r, p),y)) : -Inf    
end


