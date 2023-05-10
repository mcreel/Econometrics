# Estimates the Klein consumption equation by GMM
using Econometrics, DelimitedFiles, Statistics
function KleinMoments(θ, y, x, z)
	e = y - x*θ
	m = e.*z
end
	
function main()
cd(@__DIR__)
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
# GMM estimation using CUE
θstart = x\y # ols start values
names = ["Constant", "Profits", "Profits-1", "Wages"]
moments = θ -> KleinMoments(θ, y, x, exogs)
# estimation results using CUE
gmmtitle = "Klein model 1 GMM example CUE"
gmmresults(moments, θstart, "",  gmmtitle, names)
nothing
end
main()
