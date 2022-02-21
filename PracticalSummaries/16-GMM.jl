# Practical Summary: Ch. 16, GMM

## DGP: Poisson
using Distributions
n = 30
x = [ones(n) rand(n)]
β⁰ = [1., 2.]
λ = β -> exp.(x*β)
y = rand.(Poisson.(λ(β⁰)))

## ML, for reference
using Econometrics, SpecialFunctions
model = β -> -λ(β) + y.*x*β - loggamma.(y .+ 1.)    # note: loggamma(y+1) = log(factorial(y))
mleresults(model, β⁰);                               # but won't overflow

## GMM
using Econometrics
moments1 = β -> x.*(y - λ(β))
gmmresults(moments1, β⁰);
# note that the GMM estimator is identical to the ML estimator. You should
# be able to prove that result analytically. Because this GMM estimator is
# equivalent to ML, it is asymptotically efficient. Adding moment conditions
# won't improve asymptotic efficiency, even if they're valid.

## Trying a different GMM estimator, based on E(y)=λ, so E(y/λ)=1.
using Econometrics
moments2 = β -> x.* (y./λ(β) .- 1.0) 
gmmresults(moments2, β⁰);

## Try out overidentified GMM estimator, using both sets of moments
using Econometrics, ForwardDiff
moments3 = β -> [moments1(β) moments2(β)]
βhat, objv, V, D, W, convergence = gmmresults(moments3, β⁰);

## how to get D and Omega (though gmmresults will also give them)
avgmoments = β -> (1/n)*[x'*(y - λ(β)) x'*(y./λ(β) .- 1.)]
D = ForwardDiff.jacobian(avgmoments,βhat)'
m = moments(βhat) # the moment contributions, evaluated at estimate
Ωhat = NeweyWest(m)

## let's see how the Hansen-Sargan test can detect 
# incorrect moments
η = 0.01 # if this is different from zero, moments4, and, thus moments5, will not be valid
moments4 = β -> x.* (y./λ(β) .- 1.0 .+ η)
moments5 = β -> [moments1(β) moments4(β)]
gmmresults(moments5, β⁰);

