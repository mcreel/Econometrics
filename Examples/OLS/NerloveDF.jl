# shows some data analysis using DataFrames
using Econometrics, CSV, DataFrames, GLM
nerlove = DataFrame(CSV.File(joinpath(@__DIR__,"../Data/nerlove.csv")))
f = @formula(log(cost) ~ 1 + log(output) + log(labor) + log(fuel) + log(capital))
println()
println("using the GLM package")
@show lm(f, nerlove)
println()
println("using ols from this package")
ols(f, nerlove)
nothing
