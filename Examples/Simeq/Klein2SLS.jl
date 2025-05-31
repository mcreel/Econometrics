## Estimates the Klein consumption equation by 2SLS

## read in the data and create needed variables	
using Econometrics, CSV, DataFrames, DataFramesMeta, Statistics, Term
cd(@__DIR__)
klein = CSV.read("klein.data", DataFrame; header=true)
## construct missing lags, and drop first row that has missing data
klein.lagP = vcat(missing, [klein.P[i-1] for i=2:size(klein,1)])
klein.lagX = vcat(missing,[klein.X[i-1] for i=2:size(klein,1)])
klein.WAGES = klein.WP + klein.WG
klein = dropmissing(klein)

# define instruments
n = size(klein,1)
exogs = Matrix(@select(klein, :YEAR, :WG, :G, :T, :lagP,:lagX))
exogs = [ones(n,1) exogs]

# CONSUMPTION
println("CONSUMPTION EQUATION")
# define variables in consumption equation
y = Matrix(@select(klein, :C))
x = Matrix(@select(klein, :P, :lagP, :WAGES))
x = [ones(n,1) x]
# 2SLS estimation
names = ["Constant", "Profits", "Profits-1", "Wages"]
tsls(y, x, exogs; names=names)
printstyled("\nNote that 2SLS does not give the same results as does GMM.\n", color=:green)
printstyled("This is because the equation is overidentified. GMM is\n", color=:green)
printstyled("using an estimate of the efficient weight matrix. 2SLS does not.\n", color=:green)

