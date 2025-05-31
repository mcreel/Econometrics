# this computes the GMM estimator by SA minimization, for all of the
# Monte Carlo data sets
using Econometrics, SolveDSGE, CSV, DataFrames, Statistics, LinearAlgebra, PrettyTables
cd(@__DIR__)
include("DSGEmoments.jl")  # computes errors
include("CKlib.jl")
# needed to explore other data sets 
global const dsge = retrieve_processed_model("../GenData/CK_processed.txt")
lb, ub = PriorSupport()
θtrue = TrueParameters()
function main()
    results = zeros(1000,7)
    for rep = 1:1000
        # load the data set
        data = Matrix(CSV.read("../GenData/MCdata/mcdata-design-$rep.csv", DataFrame))
        # define CUE GMM criterion
        n = size(data,1)
        moments = theta -> sqrt(n)*DSGEmoments(theta, data)
        m = theta -> vec(mean(moments(theta),dims=1)) # 1Xg
        weight = theta -> inv(cov(moments(theta)))
        obj = theta -> m(theta)'*weight(theta)*m(theta)
        # estimate by simulated annealing
        θhat, objvalue, converged, details = samin(obj, θtrue, lb, ub; ns = 10, nt=5, verbosity = 1, rt = 0.25)
        results[rep,:] = θhat
        println("rep $rep done")
    end
    errs = results .- θtrue'
    b = mean(errs, dims=1)[:]
    m = b .+ θtrue
    s = std(errs, dims=1)[:]
    r = sqrt.(b.^2 + s.^2)[:]
    names = ["β", "γ", "ρ₁", "σ₁", "ρ₂", "σ₂", "nss"] 
    printstyled("Monte Carlo GMM for the DSGE model, 1000 reps\n", color=:green)
    pretty_table([names round.([TrueParameters() m b s r],digits=4)],
        header=["parameter", "True", "mean", "bias", "st. dev.", "rmse"]) 
end    
main()
