# try doing MCMC on GPU
# this appears to be about 40% faster than on CPU, rough guess

using PrettyTables, Pkg, CSV, Distributions, LinearAlgebra, MCMCChains, StatsPlots, CUDA
cd(@__DIR__)
Pkg.activate(".")
# defines the net and the DSGE model, and needed functions
include("Setup.jl")

function main()

net = load_trained()
Flux.testmode!(net)
net |> gpu

## Now, let's move on to Bayesian MSM using either the typical data set, or generate a new one
# load the data
data = CSV.File("dsgedata.csv") |> CSV.Tables.matrix
# transform the data the same way as was used to train net
data .-= [0.84, 0.69, 0.33, 0.05, 1.72]'
data ./= [0.51, 0.44, 0.36, 0.018, 0.34]'
X = zeros(Float32, 160, 1, 5)
X[:, 1, :] = Float32.(data)
data = tabular2conv(permutedims(Float32.(X), (3, 2, 1)))
data |> gpu
fit = net(data)
## This is the raw TCN estimate using the official data set
θnn = Float32.(UntransformParameters(fit))[:]
θnn |> gpu


#################### Define functions for MCMC ###############################

# compute mean and cov of moments, for obj fn and proposal
function simmomentscov(θ, S::Int64)
    data = Float32.(MakeData(θ, S, CKmodel))
    data |> gpu
    fit = net(data)
    m = Float32.(UntransformParameters(fit)')
    mean(m, dims=1)[:], cov(m)
end

# CUE objective, written to MAXIMIZE
@inbounds function bmsmobjective(θ, θnn, S::Int64)    # Make sure the trial parameter value is in the support
    InSupport(θ) || return -Inf
    # Compute simulated moments and covariance
    θbar, Σ = simmomentscov(θ, S)
    n = 160 # sample size
    Σ *= n * (1+1/S) # 1 for θhat, 1/S for θbar
    isposdef(Σ) || return -Inf
    err = sqrt(n)*(θnn-θbar)
    W = inv(Σ)
    Float32.(-0.5*dot(err, W, err))
end

# proposal: MVN random walk
@inbounds function proposal(current, δ, Σ)
    Float32.(rand(MvNormal(current, δ*Σ)))
end

@views function mcmc(
    θ; # TODO: prior? not needed at present, as priors are uniform
    Lₙ::Function, proposal::Function, burnin::Int=100, N::Int=1_000,
    verbosity::Int=10
)
    Lₙθ = Lₙ(θ) # Objective at data moments value
    naccept = 0 # Number of acceptance / rejections
    accept = false
    acceptance_rate = 1f0
    chain = zeros(N, size(θ, 1) + 2)
    t = time()
    for i ∈ 1:burnin+N
        θᵗ = proposal(θ) # new trial value
        Lₙθᵗ = Lₙ(θᵗ) # Objective at trial value
        # Accept / reject trial value
        accept = rand() < exp(Lₙθᵗ - Lₙθ)
        if accept
            # Replace values
            θ = θᵗ
            Lₙθ = Lₙθᵗ
            # Increment number of accepted values
            naccept += 1
        end
        # Add to chain if burnin is passed
        # @info "current log-L" Lₙθ
        if i > burnin
            chain[i-burnin,:] = vcat(θ, accept, Lₙθ)
        end
        # Report
        if verbosity > 0 && mod(i, verbosity) == 0
            tt = time()
            ttt = round(tt-t, digits=2)
            println("time for $verbosity mcmc iters: ", ttt)
            t = tt
            acceptance_rate = naccept / verbosity
            println("Current parameters (iteration i=$i)\n", round.(θ, digits=3)')
            println("acceptance_rate: ", acceptance_rate)
            naccept = 0
        end
    end
    return chain
end

#################### End Define functions for MCMC ###############################


## set up proposal and chain

# proposal
covreps = 1000
_,Σₚ = simmomentscov(θnn, covreps)
δ = 0.25 # tuning

# define objective and proposal
S = 25  # number of simulations for moments
obj = θ -> bmsmobjective(θ, θnn, S)
prop = θ -> proposal(θ, δ, Σₚ)

## run the chain
chain = mcmc(θnn, Lₙ=obj, proposal=prop, N=200)
# report results
chn = Chains(chain[:,1:end-2], ["β", "γ", "ρ₁", "σ₁", "ρ₂", "σ₂", "nss"])
plot(chn)
savefig("chain.png")
display(chn)
pretty_table([TrueParameters() θnn mean(chain[:,1:end-2],dims=1)[:]], header = (["θtrue", "θnn", "θmcmc"]))

end
main()
