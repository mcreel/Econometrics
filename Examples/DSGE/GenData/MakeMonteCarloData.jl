using SolveDSGE, Plots, DataFrames, CSV
cd(@__DIR__)
include("CKlib.jl")

if !isfile("./CK_processed.txt") 
    process_model("./CK.txt")
end    
global const dsge = retrieve_processed_model("./CK_processed.txt")

function MakeMonteCarlo()
    reps = 1000
    seed = 9999
    data = dgp(TrueParameters(), dsge, reps, seed)
    Threads.@threads for r = 1:reps
        df = DataFrame(data[r], ["y", "c", "n","r","w"])
        CSV.write("./MCdata/mcdata-design-$r.csv", df)
    end    
end

