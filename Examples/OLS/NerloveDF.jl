# shows some data analysis using DataFrames
using CSV, DataFrames, GLM
dir = dirname(dirname(pathof(Econometrics)))
nerlove = CSV.read(dir*"/Examples/Data/nerlove.csv")
f = @formula(log(cost) ~ 1 + log(output) + log(labor) + log(fuel) + log(capital))
println()
println("using the GLM package")
@show lm(f, nerlove)
println()
println("using ols from this package")
ols(f, nerlove)
nothing
