# Estimates the Klein consumption equation by GMM
using DelimitedFiles, Statistics
function KleinMoments(theta, y, x, z)
	e = y - x*theta
	m = e.*z
end
	
function main()
data = readdlm("klein.data")
# construct missing lags, and drop first row that has missing data
profits = data[:,3]
output = data[:,7]
data = [data lag(profits,1) lag(output,1)]
data = data[2:end,:]
n = size(data,1)
# define instruments
exogs = [1, 6, 8, 9, 10, 11, 12]
exogs = data[:,exogs]
exogs = [ones(n,1) exogs]
# CONSUMPTION
println("CONSUMPTION EQUATION")
# define variables in consumption equation
y = data[:,2]
profits = data[:,3]
lagprofits = data[:,11]
wp = data[:,4]
wg = data[:,8]
wages = wp + wg
# regressors in consumption equation
x = [profits lagprofits wages]
x = [ones(n,1) x]
# GMM estimation
theta = x\y
weight = 1.0
names = ["Constant", "Profits", "Profits-1", "Wages"]
moments = θ -> KleinMoments(θ, y, x, exogs)
# initial consistent estimate: only used to get moment covariance (needed for t-stats) no screen output
theta, obj_value, D, ms, convergence = gmm(moments, vec(theta), weight)
# moment covariance assuming no autocorrelation
momentcov = cov(ms)
weight = inv(momentcov)
# estimation results using efficient weight (no autocorrelation)
gmmtitle = "Klein model 1 GMM example, plain covariance"
gmmresults(moments, theta, weight, gmmtitle, names)
# moment covariance assuming autocorrelation
# note: if there really is autocorrelation,
# then lagged endogs need to be dropped as 
# instruments. This is not done here, as this
# is just meant as an example of use of NW
# covariance estimator
momentcov = NeweyWest(ms)
weight = inv(momentcov)
# estimation results using efficient weight (NW)
gmmtitle = "Klein model 1 GMM example, NW covariance"
gmmresults(moments, theta, weight, gmmtitle, names)
return
end
main()
