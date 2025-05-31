##
using Distributions, Econometrics, Plots
function main(y, tuning)
# set prior, likelihood and proposal
Prior = θ -> pdf(LogNormal(1.0,1.0), θ)
lnL = θ -> sum(logpdf.(Ref(Exponential(θ)),y)) # the log-likelihood (not averaged!)
length = 10000 # length of the chain
burnin = 1000 # drop this many, to get rid of startup bias
Proposal = θ -> rand(LogNormal(log(θ),tuning))
# get the chain, plot posterior, and descriptive stats
chain = mcmc(1.0, length, burnin, Prior, lnL, Proposal) # start value, chain length, and burnin
p = npdensity(chain[:,1]) # nonparametric plot of posterior density
plot!(p, title="posterior density: true value = 3.0") # add a title
p2 = plot(chain[end-1000:end,1])
dstats(chain)
return p, p2
end

##
n = 100 # sample size
# sample is from exponential, prior is lognormal, proposal is random walk lognormal
y = rand(Exponential(3.0),n)
# experiment with tuning, to see how acceptance changes
tuning = 0.4
p, p2 = main(y, tuning)

##
display(p)
## 
display(p2)

