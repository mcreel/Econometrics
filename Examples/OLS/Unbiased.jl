# this is a little Monte Carlo exercise that illustrates that
# the OLS estimator is unbiased when we have strong exogeneity
using Plots
reps = 1000 # number of Monte Carlo reps.
n = 20 # sample size
sig = 3.0  # st. dev. of errors

x = [ones(n,1) randn(n,1)]  # x is fixed over repeated samples
beta = [1.0, 2.0] # true beta
PopRegLine = x*beta

e = sig*randn(n,reps) # reps will be in columns
y = PopRegLine .+ e

betas = inv(x'x)*x'y

betas[2,:] .-=  2.0
histogram(betas[2,:], label = "")
gui()
#savefig("Unbiased.png")
	
