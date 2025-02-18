# This does estimation via simulated annealing, only want point estimates.

using PrettyTables, Pkg, CSV, Distributions, LinearAlgebra, CUDA

cd(@__DIR__)
Pkg.activate(".")
# defines the net and the DSGE model, and needed functions
include("Setup.jl")
include("samin.jl")

function main()

maxevals = 300
reps = 100 # number of Monte Carlo reps
# holders for results
thetahat_nn = [zeros(7) for _=1:reps] 
thetahat_sa = [zeros(7) for _=1:reps]

net = load_trained()
Flux.testmode!(net)
net |> gpu

## Now, let's move on to Bayesian MSM using either the typical data set, or generate a new one
#=
# load the data
data = CSV.File("dsgedata.csv") |> CSV.Tables.matrix
# transform the data the same way as was used to train net
data .-= [0.84, 0.69, 0.33, 0.05, 1.72]'
data ./= [0.51, 0.44, 0.36, 0.018, 0.34]'
X = zeros(Float32, 160, 1, 5)
X[:, 1, :] = Float32.(data)
data = tabular2conv(permutedims(Float32.(X), (3, 2, 1)))
=#

# compute mean and cov of moments, for obj fn and proposal
function simmoments(θ, S::Int64)
    data = Float32.(MakeData(θ, S, CKmodel))
    data |> gpu
    fit = net(data)
    m = UntransformParameters(fit)'
    mean(m, dims=1)[:]
end

# CUE objective, written to MINIMIZE
@inbounds function bmsmobjective(θ, θnn, Weight, S::Int64)    # Make sure the trial parameter value is in the support
    InSupport(θ) || return -Inf
    # Compute simulated moments and covariance
    θbar = simmoments(θ, S)
    err = θnn-θbar
    160.0 * dot(err, err)
end

# the main Monte Carlo loop
for rep = 1:reps

# Generate data
data = Float32.(MakeData(TrueParameters(), 1, CKmodel))
data |> gpu
fit = net(data)
## This is the raw TCN estimate using the official data set
θnn = UntransformParameters(fit)[:]
thetahat_nn[rep] = θnn

# get the net fit
θnn  = Float64.(θnn)
@show θnn
# get the weight matrix
data = Float32.(MakeData(θnn, 1000, CKmodel))
data |> gpu
fit = net(data)
m = UntransformParameters(fit)'
c = Float64.(cov(m))
Weight = inv(c)
@show round.(Weight, digits = 4)
# define objective
S = 100  # number of simulations for moments
obj = θ -> bmsmobjective(θ, θnn, Weight, S)
lb, ub = PriorSupport()
lb = Float64.(lb)
ub = Float64.(ub)
θhat, junk, junk, junk =  samin(obj, θnn, lb, ub, rt=0.8, nt=3, ns=1, maxevals=maxevals, functol=1e-5, paramtol=1e-3, coverage_ok=1, verbosity=3)
thetahat_sa[rep] = θhat
end

err_nn = [thetahat_nn[i] - TrueParameters() for i = 1:reps] 
err_sa = [thetahat_sa[i] - TrueParameters() for i = 1:reps] 

err_nn, err_sa, thetahat_nn, thetahat_sa
end
err_nn, err_sa, thetahat_nn, thetahat_sa = main()
mean(err_nn)
mean(err_sa)
