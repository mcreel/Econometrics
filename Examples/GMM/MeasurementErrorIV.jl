# this script does a Monte Carlo that applies IV to solve inconsistency of OLS estimator
# when there is measurement error of the regressors
using Econometrics, Plots

function GIVmoments(theta, data)
	data = [data lags(data,2)]
    data = data[3:end,:] # get rid of missings
	n = size(data,1)
	y = data[:,1]
	ylag = data[:,2]
	x = data[:,3]
	xlag = data[:,6]
	xlag2 = data[:,9]
	X = [ones(n,1) ylag x]
	e = y - X*theta
	Z = [ones(n,1) x xlag xlag2]
	m = e.*Z
end
	
function main()
# do the Monte Carlo
n = 100
sig = 1
reps = 1000
results = zeros(reps,3)
for rep = 1:reps
	x = randn(n) # an exogenous regressor
	e = randn(n) # the error term
	ystar = zeros(n)
	# generate the dep var
	for t = 2:n
	  ystar[t] = 0.0 + 0.9*ystar[t-1] + 1.0*x[t] + e[t]
	end
    # add measurement error
	y = ystar + sig*randn(n)
	# now do GMM, using the data with meas. error in both dep. var. and regressor
	ylag = lag(y,1)
	data = [y ylag x];
    data = data[2:end,:] # drop first obs, missing due to lag
	theta = [0, 0.9, 1]
	weight = 1
    # note to self: this is very slow, because it uses the general GMM method, instead of analytic GIV
    # also, the weight matrix is not optimal
    moments = theta-> GIVmoments(theta, data)
	b, junk, junk, junk, junk = gmm(moments, theta, weight)
	b = b - [0.0, 0.9, 1.0] # subtract true values, so mean should be approx. zero if consistent
    results[rep,:] = b
end

histogram(results[:,1],nbins = 30)
#savefig("givconstant.png")
histogram(results[:,2],nbins = 30)
#savefig("givrho.png")
dstats(results)	
return
end
main()

