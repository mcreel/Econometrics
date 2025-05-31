# Estimates the Klein consumption equation by GMM
using Econometrics, CSV, DataFrames, DataFramesMeta, Statistics, Term

function KleinMoments(θ, y, x, z)
	e = y - x*θ
	m = e.*z
end
	
function main()
## read in the data and create needed variables	
cd(@__DIR__)
klein = CSV.read("klein.data", DataFrame; header=true)
## construct missing lags, and drop first row that has missing data
klein.lagP = vcat(missing, [klein.P[i-1] for i=2:size(klein,1)])
klein.lagX = vcat(missing,[klein.X[i-1] for i=2:size(klein,1)])
klein.WAGES = klein.WP + klein.WG
klein = dropmissing(klein)

# CONSUMPTION
println("CONSUMPTION EQUATION")
# define instruments
n = size(klein,1)
exogs = Matrix(@select(klein, :YEAR, :WG, :G, :T, :lagP,:lagX))
exogs = [ones(n) exogs]
# define variables in consumption equation
y = Matrix(@select(klein, :C))
x = Matrix(@select(klein, :P, :lagP, :WAGES))
x = [ones(n) x]

# GMM estimation
theta = x\y  # ols start values
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
