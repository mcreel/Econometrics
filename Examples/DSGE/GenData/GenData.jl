using SolveDSGE, StatsPlots, DelimitedFiles
include("CKlib.jl")
data = CKdgp(TrueParameters(), dsge)
plot(data, legend=:outertopright, label=["output" "cons" "hours" "r" "w"])
#savefig("dsgedata.svg")
#writedlm("dsgedata.txt", data)

