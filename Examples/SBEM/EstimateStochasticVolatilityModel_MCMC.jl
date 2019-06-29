using SV, Econometrics, LinearAlgebra, Statistics, DelimitedFiles

function main()
    # load the example data
    y = readdlm("svdata.txt")
    y = y[:]
    n = size(y,1)
    m =  sqrt(n)*aux_stat(y) # generate the sample and save the data
    # these are the true params
    σe = exp(-0.736/2.0)
    ρ = 0.9
    σu = 0.363
    θtrue = [σe, ρ, σu] # true param values, on param space
    n = 1000 # sample size
    burnin = 100
    S = 100 # number of simulations
    # or generate some new data, if you prefer
    #shocks_u = randn(n+burnin,1)
    #shocks_e = randn(n+burnin,1)
    # m = SVmodel(σe, ρ, σu, n, shocks_u, shocks_e, true)
    # set up MCMC
    shocks_u = randn(n+burnin,S) # fixed shocks for simulations
    shocks_e = randn(n+burnin,S) # fixed shocks for simulations
    tuning = [0.01, 0.01, 0.01] # fix this somehow
    lb = [0.0, 0.0, 0.0]
    ub = [3.0, 0.99, 3.0]
    θinit = (ub+lb)./2.0
    lnL = θ -> logL(θ, m, n, shocks_u, shocks_e)
    Prior = θ -> prior(θ, lb, ub) # uniform, doesn't matter
    # define things for MCMC
    verbosity = true
    burnin = 100
    ChainLength = 1000
    # initial proposal moves one at a time
    Proposal = θ -> proposal1(θ, tuning, lb, ub)
    chain = mcmc(θinit, ChainLength, burnin, Prior, lnL, Proposal, verbosity)
    # keep every 10th
    i = 1:size(chain,1)
    keep = mod.(i,10.0).==0
    chain = chain[keep,:]
    # now use a MVN random walk proposal 
    Σ = cov(chain[:,1:3])
    tuning = 1.0
    for j = 1:2
        P = (cholesky(Σ)).U
        Proposal = θ -> proposal2(θ,tuning*P, lb, ub)
        θinit = vec(mean(chain[:,1:3],dims=1))
        if j == 2
            ChainLength = 10000
        end    
        chain = mcmc(θinit, ChainLength, 0, Prior, lnL, Proposal, verbosity)
        accept = mean(chain[:,end])
        #println("tuning: ", tuning, "  acceptance rate: ", accept)
        if accept > 0.35
            tuning *= 1.5
        elseif accept < 0.25
            tuning *= 0.5
        end
        # keep every 10th
        i = 1:size(chain,1)
        keep = mod.(i,10.0).==0
        θinit = vec(mean(chain[:,1:3],dims=1))
        Σ = 0.5*Σ + 0.5*cov(chain[:,1:3])
    end
    # keep every 10th to reduce autocorrelation
    i = 1:size(chain,1)
    j = mod.(i,10.0).==0
    chain = chain[j,:] 
    # plain MCMC fit
    posmean = vec(mean(chain[:,1:3],dims=1))
    inci = zeros(3)
    for i = 1:3
        lower = quantile(chain[:,i],0.05)
        upper = quantile(chain[:,i],0.95)
        inci[i] = θtrue[i] >= lower && θtrue[i] <= upper
    end
    p1 = npdensity(chain[:,1]) # example of posterior plot
    p2 = npdensity(chain[:,2]) # example of posterior plot
    p3 = npdensity(chain[:,3]) # example of posterior plot
    display(plot(p1,p2,p3))
    prettyprint([posmean inci])
    return chain
end
main()
