# Practical Review: Extremum Estimators (Ch. 13)

##
#= 	Let's start with a simple problem, where we can get the the answers using
	analytic methods. This will let us apply numeric methods, and we will be
	able to verify that they give us the correct results. To start, let's
	define some data, and do the basic OLS computations
=#	
using Econometrics
n = 5
x = [ones(n) 1:n]
β = [10., -1.]
ϵ = randn(n)
y = x*β + ϵ
ols(y,x);

##

#= 	OK, now let's do this as we would if we did not have access to
	the analytic results. First, we compute the extremum estimator
	using numeric minimization.
=#
using Statistics, Optim
objᵢ = β -> (y-x*β).^2.  	# iᵗʰ obs contrib to sum of squares
obj = β -> mean(objᵢ(β))	# average objective
tol = 1e-08
βhat = Optim.optimize(obj, zeros(2), LBFGS(), 
        	Optim.Options(g_tol = tol,x_tol=tol,f_tol=tol);
			autodiff=:forward).minimizer

##

#=	How do we get the estimated standard errors, using extremum estimation
	theory? We need estimators of I, the covariance of the score 
	contributions, and J, the Hessian matrix.
=#

# the calculations to get I hat

# first, we need the score contributions, the matrix that collects
# the derivatives of each observation's contribution to the objective

# Note to self: remember crtl-enter will evaluate a single line.

using ForwardDiff
sc =  ForwardDiff.jacobian(objᵢ, βhat) # get the score contributions
# from theory, we know that sc = -2x.*ϵcat. Let's use this verify
# that the automatic differentiation worked.
sc - (-2x.*(y - x*βhat))

##

# Next, we use the score contributions to estimate I
Ihat = zeros(2,2)
for i = 1:n
	Ihat .+= sc[i,:]*sc[i,:]'
end
Ihat ./= n
# you could also use simply cov(sc), which converges to the same thing.


##

# next, we need the estimate of the limiting Hessian, which is just
# the Hessian of the objective function, at the estimate
Jhat = ForwardDiff.hessian(obj, βhat)

# We know that this should be 2x'x/n. Let's check:
2x'x/n
##

# now, we can compute the estimated standard errors, to compare to what we saw
# from the analytic results, above
using LinearAlgebra
v∞ = inv(Jhat)*Ihat*inv(Jhat) # this is the estimate of the limiting var of √n(β-β⁰) 
se = sqrt.(diag(v∞/n))   # to get small sample est. variance, divide by n

##

# the last problem is for a correctly specified model, Case I in the notes.
# Let's look at an incorrectly specified model, Case III in the notes, to
# verify that the extremum theory works here, too.

# Let's work with Problem 1, in the exercises at the end of Chapter 13.
# In this problem, the true model is quadratic, but we erroneously estimate
# a linear model. The problem asks for an analytic solution. Here, let's
# approach it numerically.
# generate data
using Plots
function dgp(n)
	x = sort(rand(n))
	y = 1. .- x.^2 + randn(n)
	x,y
end
x,y = dgp(100)
scatter(x,y, legend=false, title="y = 1 - x² + N(0,1)")

##

# Here's an OLS fit, just to have a look
using Econometrics
n = 100
x,y = dgp(n)
X = [ones(n) x] # define regressor matrix for linear approximation about 0
b, junk = ols(y,X)
fitted = X*b
plot!(x,fitted)

##

# We know that the pseudo-true parameter values are 7/6 and -1. Let's verify
# that the OLS estimates are asymptotically normally distributed about these
# values. To do this, we will construct asymptotic 100x(1-α)% confidence intervals,
# and verify that the pseudo-true values lie inside them approximately 100x(1-α)% of
# the times we repeat the procedure, at least when n is large enough
using LinearAlgebra, Statistics, Distributions
n = 20 				# try small and large values here,
					# to see accuracy of asymptotic approximation
reps = 10000
α = 0.05			# try out 90, 95 and 99% CIs
crit = quantile(Normal(),1-α/2)
β⁰ = [7/6; -1.]
inci = zeros(reps, 2)
for i = 1:reps
	x,y = dgp(n)
	X = [ones(n) x] # define regressor matrix for linear approximation about 0
	βhat, vβhat, junk = ols(y, X, silent=true) # note: OLS uses the sandwich covariance estimator, by default
	se = sqrt.(diag(vβhat))
	inci[i,:] = (β⁰ .>= βhat .- crit*se) .& (β⁰ .<= βhat .+ crit*se)
end
ci = Int64(100*(1-α))
println("Coverage of $ci% CIs")
mean(inci, dims=1) # these should be approximately 1-α, at least when n is large enough

##

# To conclude, this summary has shown that the extremum estimation theory is
# working for two cases, a correctly specified model, and an incorrectly specified
# model. We have seen how to get I and J, now to compute the asymptotic variance,
# and how to verify that confidence intervals (which are computed using I and J)
# have correct coverage, at least asymptotically.


