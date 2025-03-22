# This does estimation via simulated annealing, only want point estimates.

using Pkg, CSV, DataFrames, Statistics, PrettyTables, Term
include("Setup.jl")

cd(@__DIR__)
Pkg.activate(".")

#function main()
infile = "PointEstimationResults.txt"

ind = [2,5]
θtrue = TrueParameters()[ind]
data = Matrix(CSV.read(infile, DataFrame, header=false))

θnn = data[:,ind]
θsa = data[:, ind .+ 7]

err_nn = (θnn .- θtrue')
err_sa = (θsa .- θtrue')

b_nn =  mean(err_nn, dims=1)'
b_sa = mean(err_sa, dims=1)'
r = size(data,1)
println("\n____________________________________")
println(@green "Results based on $r replications")
println("Bias")
names = ["γ", "ρ_η"]   
pretty_table([names b_nn b_sa], header=["param.", "θnn", "θsa"])

rmse_nn =  sqrt.(mean(err_nn .^2, dims=1)')
rmse_sa = sqrt.(mean(err_sa .^2, dims=1)')
println("RMSE")
pretty_table([names rmse_nn rmse_sa], header=["param", "θnn", "θsa"])



#end
#b = main()
