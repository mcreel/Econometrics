# This does MLE and then MCMC, either using raw statistic, or using NN transform,
# depending on the argument usenn
using Flux, Econometrics, LinearAlgebra, Statistics, DelimitedFiles
using BSON:@load
include("transform.jl")
include("lnL.jl")

# uniform random walk in one dimension
function proposal1(current, tuning)
    trial = copy(current)
    i = rand(1:size(trial,1))
    trial[i] += tuning[i].*randn()
    return trial
end

# MVN random walk, or occasional draw from prior
function proposal2(current, cholV)
    current + cholV'*randn(size(current))
end

function MCMC(m, usenn, info, useJacobian=true)
    S = 100 # number of simulations
    lb, ub = PriorSupport()
    nParams = size(lb,1)
    # get the trained net
    if usenn
        @load "best.bson" model
        m = transform(m', info)
        m = Float64.(model(m'))
        θinit = PriorMean() # prior mean as initial θ
        lnL = θ -> LL(θ, m, S, model, info, useJacobian)
    else          
        θinit = PriorMean() # prior mean as initial θ
        lnL = θ -> LL(θ, m, S, useJacobian)
    end
    # use a rapid SAMIN to get good initialization values for chain
    obj = θ -> -1.0*lnL(θ)
    θmile, junk, junk, junk = samin(obj, θinit, lb, ub; coverage_ok=0, maxevals=100000, verbosity = 0, rt = 0.5)
    prior = θ -> Prior(θ) # uniform, doesn't matter
    # define things for MCMC
    verbosity = false
    ChainLength = 1000
    MCMCburnin = 100
    tuning = 0.2/sqrt(12.0)*(ub-lb) # two tenths of a standard. dev. of prior
    Proposal = θ -> proposal1(θ, tuning)
    chain = mcmc(θmile, ChainLength, MCMCburnin, prior, lnL, Proposal, verbosity)
    # now use a MVN random walk proposal with updates of covariance and longer chain
    # on final loop
    Σ = NeweyWest(chain[:,1:nParams])
    tuning = 1.0
    MC_loops = 5
    @inbounds for j = 1:MC_loops
        P = try
            P = (cholesky(Σ)).U
        catch
            P = diagm(diag(Σ))
        end    
        Proposal = θ -> proposal2(θ,tuning*P)
        if j == MC_loops
            ChainLength = 10000
        end    
        θinit = mean(chain[:,1:nParams],dims=1)[:]
        chain = mcmc(θinit, ChainLength, 0, prior, lnL, Proposal, verbosity)
        if j < MC_loops
            accept = mean(chain[:,end])
            if accept > 0.35
                tuning *= 1.5
            elseif accept < 0.25
                tuning *= 0.25
            end
            Σ = 0.5*Σ + 0.5*NeweyWest(chain[:,1:nParams])
        end    
    end
    chain = chain[:,1:nParams]
    return chain, θmile
end
