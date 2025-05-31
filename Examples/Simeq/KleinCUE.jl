# Estimates the Klein consumption equation by CUE GMM
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
