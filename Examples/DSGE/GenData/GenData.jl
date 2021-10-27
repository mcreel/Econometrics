using SolveDSGE, Plots, DelimitedFiles, DataFrames, CSV
include("CKlib.jl")
data = dgp(TrueParameters(), dsge, 1)[1]
plot(data, legend=:outertopright, label=["output" "cons" "hours" "r" "w"])
#savefig("dsgedata.svg")
#writedlm("dsgedata.txt", data)
df = DataFrame(data, ["output", "cons", "hour","r","w"])
#CSVwrite("dsgedata.csv", df)


