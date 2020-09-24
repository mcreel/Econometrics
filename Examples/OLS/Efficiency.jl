# this is a little Monte Carlo exercise that illustrates that
# the OLS estimator is efficient. It compares the OLS estimator
# to an estimator that is the average of the OLS estimators using 
# 3 subsamples
using Plots
function main()
pyplot()    
reps = 10000 # number of Monte Carlo reps.
n = 21 # sample size
sig = 2.0  # st. dev. of errors

x = [ones(n,1) randn(n,1)]  # x is fixed over repeated samples
beta = [1.0, 2.0] # true beta
PopRegLine = x*beta

betas = zeros(reps,2) # holder for results

for i = 1:reps
	# generate dep var
	e = sig*randn(n)
	y = PopRegLine + e
	# the OLS estimator
	betas[i,1] = (y\x)[2]
	# the average of split sample estimator
	y1 = y[1:7,:]
	y2 = y[8:14,:]
	y3 = y[15:21,:]
	x1 = x[1:7,:]
	x2 = x[8:14,:]
	x3 = x[15:21,:]
	beta_ss1 = y1\x1
	beta_ss2 = y2\x2
	beta_ss3 = y3\x3
	beta_ss = (beta_ss1 + beta_ss3 + beta_ss3)/3.0 # average the 3 estimators
	betas[i,2] = beta_ss[2]
end

histogram(betas[:,1], label="", nbins=30, show=true)
#savefig("efficiency1.svg")
histogram(betas[:,2], label="", nbins=30, reuse=false, show=true)
#savefig("efficiency2.svg")
return
end
main()
	
