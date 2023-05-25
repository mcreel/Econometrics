##
using Econometrics, LinearAlgebra, Statistics, DelimitedFiles, SolveDSGE, MCMCChains, Plots, Distributions
include("DSGEmoments.jl")
include("DSGEmodel.jl") # defines prior and log-likelihood
include("CKlib.jl")

function main()
    cd(@__DIR__)
    data = readdlm("dsgedata.txt")
    lb, ub = PriorSupport()
    θtrue = TrueParameters()
    # start values from GMM
    θinit = [
    0.9900485808486286,
    2.0271722237817826,
    0.9013513634634782,
    0.01901778177965745,
    0.7131947676648516,
    0.010190093131189638,
    0.333600599933210]
    lnL = θ -> logL(θ, data)
    # define things for MCMC
    verbosity = true
    burnin = 10000
    ChainLength = 10000
    # Initial chain to get covariance of proposal
    Proposal = θ -> proposal1(θ, tuning)
    tuning = [0.0001, 0.10, 0.05, 0.005, 0.005, 0.005, 0.01]
    chain = mcmc(θinit, ChainLength, burnin, Prior, lnL, Proposal, verbosity)
    # now use a MVN random walk proposal, and adjust tuning in a loop 
    Σ = cov(chain[:,1:7])
    tuning = 0.05
    for j = 1:5
        Proposal = θ -> rand(MvNormal(θ,tuning*Σ))
        if j == 5
            ChainLength = 100000
        end    
        chain = mcmc(θinit, ChainLength, burnin, Prior, lnL, Proposal, verbosity)
        accept = mean(chain[:,end])
        if accept > 0.3
            tuning *= 1.5
        elseif accept < 0.2
            tuning /= 1.5
        end
        θinit = vec(mean(chain[:,1:7],dims=1))
        Σ = NeweyWest(chain[:,1:7])
    end
    # plain MCMC fit
    posmean = vec(mean(chain[:,1:7],dims=1))
    inci = zeros(7)
    for i = 1:7
        lower = quantile(chain[:,i],0.05)
        upper = quantile(chain[:,i],0.95)
        inci[i] = θtrue[i] >= lower && θtrue[i] <= upper
    end
    prettyprint([posmean inci])
    return chain
end

##
chain = main()

## visualize results
p = npdensity(chain[:,2]) # example of posterior plot
display(p)
savefig("gamma.svg")

##
chn = Chains(chain[:,1:7], ["β", "γ", "ρ₁", "σ₁", "ρ₂", "σ₂", "nss"])
display(chn)
display(plot(chn))
savefig("allparams.svg")


