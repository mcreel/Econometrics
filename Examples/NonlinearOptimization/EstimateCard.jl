# exercise estimating a linear model by iterative minimization,
# and verifying that OLS gives same results.
# If the data is missing, run BasicDataAnalysis.jl
# in the Examples/Julia directory

using CSV, DataFrames, StatsModels, Econometrics
# prepare the data
card = CSV.read("../Julia/cooked.csv", DataFrame)
display(card)

y = card[:,1]
x = Matrix{Float64}([ones(size(card,1)) card[:,2:7]])

# define the objective function and start value
obj = theta -> (y-x*theta)'*(y-x*theta)
startval = zeros(size(x,2))
# do the minimization
thetahat, objvalue = fminunc(obj, startval) 
println("the OLS estimates by numeric min: obj. value: ", round(objvalue,digits=5))
prettyprint(thetahat)

# verify by using OLS, which uses the analytic solution
println("verifying by OLS:")
ols( @formula(lnwage ~ 1 + educ + exper + expsq + black + south + smsa), card)
nothing
