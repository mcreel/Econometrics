# use this to time regular mcmc(), but without plots.
using Distributions, Random
function main()
    y = rand(Exponential(3.0),30)
    # set prior, likelihood and proposal
    Prior = θ -> pdf.(Ref(LogNormal(1.0,1.0)), θ)
    lnL = θ -> sum(logpdf.(Ref(Exponential(θ)), y))
    tuning = 0.5
    Proposal = θ -> rand(LogNormal(log(θ),tuning))
    # get the chain, plot posterior, and descriptive stats
    chain = mcmc(1.0, 100000, 10000, Prior, lnL, Proposal, false) # start value, chain length, and burnin 
end


