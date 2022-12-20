using SolveDSGE, Plots, DelimitedFiles, DataFrames, CSV
cd(@__DIR__)
include("CKlib.jl")

if !isfile("./CK_processed.txt") 
    process_model("./CK.txt")
end    
global const dsge = retrieve_processed_model("./CK_processed.txt")

function GenData()
# this block reads and processes the file, leave it be
data = dgp(TrueParameters(), dsge, 1)[1]
display(plot(data, legend=:outertopright, label=["output" "cons" "hours" "r" "w"]))
#savefig("dsgedata.svg")
#writedlm("dsgedata.txt", data)
df = DataFrame(data, ["output", "cons", "hour","r","w"])
display(df)
#CSVwrite("dsgedata.csv", df)
return true
end

