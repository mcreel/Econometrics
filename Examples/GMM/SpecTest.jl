# this script does a Monte Carlo that applies IV to solve inconsistency of OLS estimator
# when there is measurement error of the regressors
using Econometrics, Distributions, LinearAlgebra, Random

function GIVmoments(theta, y, x, inst)
	e = y - x*theta
	ms = e.*inst
end

function main()
# do the Monte Carlo
n = 100
sig = 1
reps = 1000
results = zeros(reps,3)
ystar = zeros(n)
for rep = 1:reps
	x = randn(n) # an exogenous regressor
	e = randn(n) # the error term
	# generate the dep var
	for t = 2:n
	  ystar[t] = 0.0 + 0.9*ystar[t-1] + 1.0*x[t] + e[t]
	end
    # add measurement error
	y = ystar + sig*randn(n)
	# now do GMM, using the data with meas. error in both dep. var. and regressor
	data = [y lag(y,1) x lags(x,2)];
    data = data[3:end,:] # drop first obsns, missing due to lags
    y = data[:,1]
    ylag = data[:,2]
    x = data[:,3]
    xlags = data[:,4:end]
    regressors = [ones(size(y)) ylag x]
	instruments = [ones(size(y)) x xlags]
    xhat = instruments*(instruments\regressors)
    thetahat = xhat\y
    ms = GIVmoments(thetahat, y, regressors, instruments)
    weight = inv(cov(ms))
    moments = theta-> GIVmoments(theta, y, regressors, instruments)
	junk, objvalue, junk, junk, junk = gmm(moments, thetahat, weight)
   	objvalue = size(y,1)*objvalue
	df = size(ms,2) - size(thetahat,1)
	reject = objvalue .> quantile.(Ref(Chisq(df)),[0.9 0.95 0.99]) # does the test reject at these signif levels?
    results[rep,:] = reject
end
cnames = ["10%", "5%", "1%"]
println("rejection frequencies, nominal 10% 5% and 1%:")
prettyprint(mean(results,dims=1),cnames)
return
end
main()
