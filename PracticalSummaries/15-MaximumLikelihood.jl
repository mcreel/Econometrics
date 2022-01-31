# Practical Summary: Ch. 15, Maximum likelihood 

## Draw some data from the exponential distribution
using Distributions, Econometrics, Plots
θ⁰ = 3.0
n = 100
y = rand(Exponential(θ⁰),n)
histogram(y, normalize=true)
plot!(Exponential(θ⁰))


## Here are ML resutls, for reference
using Econometrics, Distributions
logℒᵢ = θ -> logpdf.(Exponential(θ[1]),y)
θstart = [1.]
θhat = mleresults(logℒᵢ, θstart, y)

##
# for this model, you should be able to work out, analytically,
# that the ML estimator is simply the sample mean. To confirm that:
mean(y)

## Let's get the results from basic theory: first the maximizer of log likelihood
using Optim
tol = 1e-08
s = θ -> -mean(logℒᵢ(θ))  # negative average logℒ
θhat = Optim.optimize(s, θstart, LBFGS(), 
                            Optim.Options(
                            g_tol = tol,
                            x_tol=tol,
                            f_tol=tol);autodiff=:forward).minimizer

## Now get the t-stats, from ℐhat and 𝒥  hat.
using ForwardDiff
sc =  ForwardDiff.jacobian(logℒᵢ, θhat) # get the score contributions
ℐhat = mean(sc.^2.)
𝒥hat = -ForwardDiff.hessian(s, θhat)
# three forms of estimated variance
V1 = inv(ℐhat)/n
se1 = sqrt(V1)
V2 = inv(-𝒥hat)/n
se2 = sqrt(V2)
V3 =  inv(𝒥hat)*ℐhat*inv(𝒥hat)/n
se3 = sqrt(V3)
##


# efficiency: compare to regression
    #constant only
    # with regressors (nonlinear least squares)    

# information matrix equality

# likelihood ratio test