## MeasurementError

## Author: Michael Creel <michael@yosemite>
## Created: 2010-03-04
## adapted to Julia 27 Aug. 2017

# this script does a Monte Carlo that shows the inconsistency of OLS estimator
# when there is measurement error of the regressors

using Plots, LinearAlgebra, Econometrics
function main()
# the function to be Monte Carlo'ed
function wrapper(n, sig)
	x = randn(n) # an exogenous regressor
	e = randn(n) # the error term
	ystar = zeros(n)
	# generate the dep var
	for t = 2:n
	  ystar[t,:] = 0.0 + 0.9*ystar[t-1] + 1.0*x[t] + e[t]
	end
    # add measurement error
	y = ystar + sig*randn(n)
	# now do OLS, using the data with meas. error in both dep. var. and regressor
	ylag = lag(y,1)
	data = [y ylag x]
	data = data[2:n,:] # drop first obs, missing due to lag
	y = data[:,1]
	ylag = data[:,2]
	x = data[:,3]
	x = [ones(size(x,1),1) ylag x] # data matrix for ols
	b = x\y
	b = b' - [0.0 0.9 1.0] # subtract true values, so mean should be approx. zero if consistent
end

# do the Monte Carlo
n = 100
sig = 0.0
reps = 1000
bs = zeros(reps,3)
for rep = 1:reps
    bs[rep,:] = wrapper(n, sig)
end
# analyze results
histogram(bs[:,2])
sig != 0.0 ? savefig("ylag_n100.svg") : savefig("ylag_n100_no_error.svg")
gui()
dstats(bs);
return
end
main()
