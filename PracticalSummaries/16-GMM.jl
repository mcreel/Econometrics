# Practical Summary: Ch. 16, GMM

## DGP: Poisson
using Distributions
n = 30
x = [ones(n) rand(n)]
β⁰ = [1., 2.]
λ = β -> exp.(x*β)  # the usual Poisson parameterization
y = rand.(Poisson.(λ(β⁰))) # draw the data
# data for problem set 3
# x = hcat(ones(10), [-1,-1,1,0,-1,-1,1,1,1,2])
# y = [0,0,0,1,1,5,8,16,20,30]
# β⁰ = [1.3, 1.0]  # use decent start values for the problem set data

## ML, for reference
using Econometrics, SpecialFunctions
# the log-likelihood
model = β -> -λ(β) + y.*x*β - loggamma.(y .+ 1.)    # note: loggamma(y+1) = log(factorial(y))
mleresults(model, β⁰);                               # but won't overflow

## GMM
using Econometrics
moments1 = β -> x.*(y - λ(β))
βhat, objv, V, D, W, convergence = gmmresults(moments1, β⁰);
# note that the GMM estimator is identical to the ML estimator. You should
# be able to prove that result analytically, by showing that the moment conditions
# are the same as the ML score vector. Because this GMM estimator is
# equivalent to ML, it is asymptotically efficient. Adding moment conditions
# won't improve asymptotic efficiency, even if they're valid.

## Trying a different GMM estimator, based on E(y)=λ, so E(y/λ)=1.
using Econometrics
moments2 = β -> x.* (y./λ(β) .- 1.0) 
βhat, objv, V, D, W, convergence = gmmresults(moments2, β⁰);


## Try out overidentified GMM estimator, using both sets of moments

## first two step
using Econometrics, ForwardDiff
moments3 = β -> [moments1(β) moments2(β)] # stacking them horizontally, so there are 4 moment conditions
W = eye(4)
βhat, objv, V, D, W, convergence = gmmresults(moments3, β⁰, W, "first step", "", false);
W = inv(cov(moments3(βhat)))
βhat, objv, V, D, W, convergence = gmmresults(moments3, β⁰, W, "second step efficient");

## now CUE for the overidentified version
using Econometrics, ForwardDiff
moments3 = β -> [moments1(β) moments2(β)]
βhat, objv, V, D, W, convergence = gmmresults(moments3, β⁰);


## how to get D and Omega (though gmmresults will also give them)
avgmoments = β -> (1/n)*[x'*(y - λ(β)) x'*(y./λ(β) .- 1.)]
D = ForwardDiff.jacobian(avgmoments,βhat)'
m = moments3(βhat) # the moment contributions, evaluated at estimate
Ωhat = NeweyWest(m)

## let's see how the Hansen-Sargan test can detect 
# incorrect moments
η = 0.02 # if this is different from zero, moments4 will not be valid
moments4 = β -> [moments1(β) .+ η moments2(β)]
βhat, objv, V, D, W, convergence = gmmresults(moments4, β⁰);

