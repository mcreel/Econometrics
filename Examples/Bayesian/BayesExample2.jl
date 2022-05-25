##
using Distributions, Econometrics, Plots
function main()
# sample is from exponential, prior is lognormal, proposal is random walk lognormal
y = rand(Exponential(3.0),30)
# set prior, likelihood and proposal
Prior = θ -> pdf(LogNormal(1.0,1.0), θ)
lnL = θ -> sum(logpdf.(Ref(Exponential(θ)),y)) # the log-likelihood (not averaged!)
tuning = 0.5
Proposal = θ -> rand(LogNormal(log(θ),tuning))
# get the chain, plot posterior, and descriptive stats
chain = mcmc(1.0, 100000, 10000, Prior, lnL, Proposal) # start value, chain length, and burnin
p = npdensity(chain[:,1]) # nonparametric plot of posterior density
plot!(p, title="posterior density: true value = 3.0") # add a title
p2 = plot(chain[end-1000:end,1])
dstats(chain)
return p, p2
end
p, p2 = main()

##
display(p)
## 
display(p2)

