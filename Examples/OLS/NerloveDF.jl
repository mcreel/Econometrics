## shows some data analysis using DataFrames
using Econometrics, CSV, DataFrames, GLM
cd(@__DIR__)
nerlove = DataFrame(CSV.File("../Data/nerlove.csv"))
f = @formula(log(cost) ~ 1 + log(output) + log(labor) + log(fuel) + log(capital))

##
println("using the GLM package")
@show lm(f, nerlove)

##
println("using ols from this package")
println("note: this uses White's st. errors, GLM uses plain vanilla")
ols(f, nerlove)
nothing
