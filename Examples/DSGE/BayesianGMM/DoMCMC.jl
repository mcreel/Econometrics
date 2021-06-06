using LinearAlgebra, Statistics, DelimitedFiles, SolveDSGE, MCMCChains, Plots
include("DSGEmoments.jl")
include("DSGEmodel.jl") # defines prior and log-likelihood
include("CKlib.jl")

function main()
    data = readdlm("dsgedata.txt")
    lb, ub = PriorSupport()
    θtrue = TrueParameters()
    θinit = (ub + lb) / 2.0 # prior mean to start
    lnL = θ -> logL(θ, data)
    # define things for MCMC
    verbosity = true
    burnin = 10000
    ChainLength = 50000
    # initial proposal moves one at a time
    Proposal = θ -> proposal1(θ, tuning)
    tuning = [0.0001, 0.1, 0.01, 0.001, 0.05, 0.001, 0.001]
    chain = mcmc(θinit, ChainLength, burnin, Prior, lnL, Proposal, verbosity)
    # now use a MVN random walk proposal 
    Σ = NeweyWest(chain[:,1:7])
    tuning = 0.1
    for j = 1:5
        P = (cholesky(Σ)).U
        Proposal = θ -> proposal2(θ,tuning*P)
        θinit = vec(mean(chain[:,1:7],dims=1))
        if j == 10
            ChainLength = 100000
        end    
        chain = mcmc(θinit, ChainLength, 0, Prior, lnL, Proposal, verbosity)
        accept = mean(chain[:,end])
        #println("tuning: ", tuning, "  acceptance rate: ", accept)
        if accept > 0.35
            tuning *= 1.5
        elseif accept < 0.25
            tuning *= 0.5
        end
        θinit = vec(mean(chain[:,1:7],dims=1))
        Σ = 0.5*Σ + 0.5*NeweyWest(chain[:,1:7])
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
chain = main()
# visualize results
p = npdensity(chain[:,2]) # example of posterior plot
display(p)
#savefig("gamma.svg")
chn = Chains(chain[:,1:7], ["β", "γ", "ρ₁", "σ₁", "ρ₂", "σ₂", "nss"])
display(chn)
#savefig("allparams.svg")


