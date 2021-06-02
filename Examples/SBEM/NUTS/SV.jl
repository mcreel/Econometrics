# simple example of method of simulated moments estimated by
# NUTS MCMC

# load packages
using DelimitedFiles, TransformVariables, LogDensityProblems, DynamicHMC, MCMCDiagnostics,
    Parameters, Statistics, Distributions, ForwardDiff, LinearAlgebra

include("SVlib.jl")

# Define a structure for the problem
# Should hold the data and  the parameters of prior distributions.
struct MSM_Problem{Tm <: Vector{Float64}, Tn <: Int, Tshocks_u <: Array{Float64,2}, Tshocks_e <: Array{Float64,2}}
    "statistic"
    m::Tm
    "sample size"
    n::Tn
    "shocks"
    shocks_u::Tshocks_u
    shocks_e::Tshocks_e

end

# Make the type callable with the parameters *as a single argument*.
function (problem::MSM_Problem)(θ)
    @unpack m, n, shocks_u, shocks_e = problem   # extract the data
    @unpack σe, ρ, σu = θ         # extract parameters (only one here)
    S = size(shocks_u,2)
    k = size(m,1)
    ms = zeros(eltype(SVmodel(σe, ρ, σu, n, shocks_u[:,1], shocks_e[:,1])), S, k)
    Threads.@threads for s = 1:S
        ms[s,:] = SVmodel(σe, ρ, σu, n, shocks_u[:,s], shocks_e[:,s])
    end
    mbar = mean(ms,dims=1)[:]
    Σ = cov(ms)
    x = (m .- mbar)
    logL = try
        logL = -0.5*log(det(Σ)) - 0.5*x'*inv(Σ)*x
    catch
        logL = -Inf
    end    
end
function main()
# generate data
σe = exp(-0.736/2.0)
ρ = 0.9
σu = 0.363
n = 1000 # sample size
burnin = 100
S = 100 # number of simulations
shocks_u = randn(n+burnin,1)
shocks_e = randn(n+burnin,1)
# load the example data
y = readdlm("../svdata.txt")
m =  sqrt(n)*aux_stat(y) # statistic
m = m[:]
# or use a new draw
#m = SVmodel(σe, ρ, σu, n, shocks_u, shocks_e)
shocks_u = randn(n+burnin,S) # fixed shocks for simulations
shocks_e = randn(n+burnin,S) # fixed shocks for simulations
# original problem, without transformation of parameters
p = MSM_Problem(m, n, shocks_u, shocks_e)
# define the transformation of parameters (in this case, priors are uniform on segments of real line)
# σu ~ U(0,3), ρ ~U(0,1), σe ~ U(0,3)
problem_transformation(p::MSM_Problem) =
        as((σe=as(Real, 0.0, 3.0), ρ=as(Real, 0.0, 1.0) ,σu=as(Real,0.0,3.0)))
# Wrap the problem with the transformation
t = problem_transformation(p)
P = TransformedLogDensity(t, p)
# use AD for the gradient
∇P = ADgradient(:ForwardDiff, P)
# Sample from the posterior. `chain` holds the chain (positions and
# diagnostic information), while the second returned value is the tuned sampler
# which would allow continuation of sampling.
n = dimension(problem_transformation(p))
chain, NUTS_tuned = NUTS_init_tune_mcmc(∇P, 1000)
# We use the transformation to obtain the posterior from the chain.
posterior = transform.(Ref(t), get_position.(chain));
return posterior, chain
end
