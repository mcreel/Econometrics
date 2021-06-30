using SolveDSGE, StatsPlots, DelimitedFiles
include("CKlib.jl")
data = dgp(TrueParameters(), dsge, 1)[1]
plot(data, legend=:outertopright, label=["output" "cons" "hours" "r" "w"])
#savefig("dsgedata.svg")
#writedlm("dsgedata.txt", data)

