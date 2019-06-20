using StatsPlots, DelimitedFiles

# generate the dats
#=
include("../SBEM/SVlib.jl")
σe = exp(-0.736/2.0)
ρ = 0.9
σu = 0.363
θtrue = [σe, ρ, σu] # true param values, on param space
n = 1000 # sample size
burnin = 100
shocks_u = randn(n+burnin,1)
shocks_e = randn(n+burnin,1)
m = SVmodel(σu, ρ, σe, n, shocks_u, shocks_e, true)
=#
y = readdlm("svdata.txt")
p1 = plot(y)
p2 = density(y)
plot(p1, p2, layout=(2,1))
gui()
savefig("svdata.svg")
dstats(y)
