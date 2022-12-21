using SolveDSGE, Plots, DelimitedFiles, DataFrames, CSV
cd(@__DIR__)
include("CKlib.jl")

if !isfile("./CK_processed.txt") 
    process_model("./CK.txt")
end    
global const dsge = retrieve_processed_model("./CK_processed.txt")

function GenData(silent=false)

#this block reads and processes the file, leave it be
data = dgp(TrueParameters(), dsge, 1)[1]
df = DataFrame(data, ["output", "cons", "hour","r","w"])
if !silent
    display(plot(data, legend=:outertopright, label=["output" "cons" "hours" "r" "w"]))
    #savefig("dsgedata.svg")
    display(df)
end    
#writedlm("dsgedata.txt", data)
#CSVwrite("dsgedata.csv", df)
return true
end

