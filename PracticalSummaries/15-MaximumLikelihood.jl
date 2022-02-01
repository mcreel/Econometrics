# Practical Summary: Ch. 15, Maximum likelihood

# First, let's work with a simple one parameter model,
# for which we can get analytic results.

## Draw some data from the exponential distribution
using Distributions, Econometrics, Plots
θ⁰ = 4.0
n = 100
y = rand(Exponential(θ⁰),n)
histogram(y, normalize=true)
plot!(Exponential(θ⁰))

## Here are the ML results, for reference
using Econometrics, Distributions
logℒᵢ = θ -> logpdf.(Exponential(θ[1]),y)
θstart = [1.]
θhat = mleresults(logℒᵢ, θstart, y)[1]

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
# when you set n large (above), you should see that ℐhat ⩬ -𝒥hat  
# also, note that 𝒥hat = 1/θhat², which is a result we can get analytically

# three forms of estimated variance
V1 = inv(ℐhat)/n
se1 = sqrt(V1)
V2 = inv(-𝒥hat)/n
se2 = sqrt(V2)
V3 =  inv(𝒥hat)*ℐhat*inv(𝒥hat)/n
se3 = sqrt(V3)
[se1 se2 se3]   # note that the estimators are a little different from one another
                # the last one, sandwich, is what's reported in mleresults.
##

# now let's work with a model that has regressors
using Distributions, Econometrics, Plots
n = 1000
x = [ones(n) randn(n) rand(n)]
β⁰ = [-0.5, 1., 1.]
θ = β ->  exp.(x*β)
y = rand.(Exponential.(θ(β⁰))) # independent, non-identically distributed observations

## Here's the exponential (correct) model
using Econometrics, Distributions, ForwardDiff, Calculus
logℒᵢ = β -> [logpdf(Exponential(exp(x[i,:]'*β)),y[i]) for i = 1:n]
s(b) = mean(logℒᵢ(b))
βstart = zeros(3)
βhat, objval, junk = mleresults(logℒᵢ, βstart, y)
sc =  ForwardDiff.jacobian(logℒᵢ, βhat) # get the score contributions
ℐhat = zeros(3,3)
for i = 1:n
	ℐhat .+= sc[i,:]*sc[i,:]'
end
ℐhat ./= n
𝒥hat = Calculus.hessian(s, βhat, :central)
ℐhat+𝒥hat # pretty close to zeros, as information matrix equality tells us

## Likelihood ratio test
using Distributions, Econometrics
# a true restriction
R = [1 1 1]
r = 1.5
obj(b) = -s(b) # need to minimize
βhatR, objvalR, junk = fmincon(obj, βhat, R, r)
objvalR = -objvalR # back to maximization
LR = 2*n*(objval - objvalR)
pval = 1. - cdf(Chisq(1), LR)
println("LR test: $R*β=$r   LR stat: $LR. p-value: $pval")
# a false restriction
R = [0 1 -1]
r = 1.
obj(b) = -s(b)
βhatR, objvalR, junk = fmincon(obj, βhat, R, r)
objvalR = -objvalR # back to maximization
LR = 2*n*(objval - objvalR)
pval = 1. - cdf(Chisq(1), LR)
println("LR test: $R*β=$r   LR stat: $LR. p-value: $pval")

## What about trying a χ² model? If we don't know the true density, we don't know which to use
logℒᵢ = β -> [logpdf(Chisq(exp(x[i,:]'*β)),y[i]) for i = 1:n]
βstart = zeros(3)
logℒᵢ(zeros(3))
βhat = mleresults(logℒᵢ, βstart, y)[1]
sc =  ForwardDiff.jacobian(logℒᵢ, βhat) # get the score contributions
ℐhat = zeros(3,3)
for i = 1:n
	ℐhat .+= sc[i,:]*sc[i,:]'
end
ℐhat ./= n
𝒥hat = Calculus.hessian(s, βhat, :central)
ℐhat+𝒥hat # a lot farther away from zeros, at least when n is large.

# How can we tell which model is best? The information criteria
# favor the exponential model, clearly. Comparing I and J,
# the information matrix equality holds much better for the
# exponential model than for the Chi-square. The IM test formalizes
# this, but here, we can see the intuitive basis for the test.


