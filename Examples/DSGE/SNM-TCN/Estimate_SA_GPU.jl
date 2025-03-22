# This does estimation via simulated annealing, only want point estimates to
# check bias correction
# Also, focussing only on gamma and rho_eta, as these are the only two with 
# significant bias. This is to speed it up.

using Pkg, CSV, Distributions, LinearAlgebra, CUDA, Tables

cd(@__DIR__)
Pkg.activate(".")
# defines the net and the DSGE model, and needed functions
include("Setup.jl")
include("samin.jl")

function main()
outfile = "PointEstimationResults.txt"


maxevals = 300
reps = 6 # number of Monte Carlo reps

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

# get the net fit
θnn  = Float64.(θnn)
@show θnn
# get the weight matrix
data = Float32.(MakeData(θnn, 1000, CKmodel))
data |> gpu
fit = net(data)
m = Float64.(UntransformParameters(fit)')
c = diagm(diag(cov(m)))
Weight = inv(c)
@show round.(Weight, digits = 4)
# define objective
S = 100  # number of simulations for moments
obj = θ -> bmsmobjective(θ, θnn, Weight, S)
ll, uu = PriorSupport()
lb = copy(θnn)
lb[2] = ll[2]
lb[5] = ll[5]
ub = copy(θnn)
ub[2] = uu[2]
ub[5] = uu[5]

θsa, junk, junk, junk =  samin(obj, θnn, lb, ub, rt=0.5, nt=3, ns=1, maxevals=maxevals, functol=1e-5, paramtol=1e-3, coverage_ok=1, verbosity=3)
a = zeros(1,14)
a[1,1:7] = θnn 
a[1,8:14] = θsa 
CSV.write(outfile, Tables.table(round.(a, digits=5)), writeheader=false, append=true) 
end
end
main()
