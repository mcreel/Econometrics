#= 
learn a bit more about MCMC. This explores using
a proposal density that takes into account the 
correlations of the posterior, to get better
mixing

* play with scaling to see effects on acceptance rate
* examine chains with different acceptance rates
* look at scatter plots to see posterior correlations
=#
using Plots, LinearAlgebra, Statistics, Econometrics
include("QIVmodel.jl")
function main()
LNW, X, Z = getdata()
n = size(LNW,1)
mcmcreps = 100000
burnin = 100000
# use st. errs. from ordinary QR as quide to setting tuning.
basetuning = [0.1, 0.05, 0.03, 0.005, 0.3, 0.2, 0.2]
# adjust this up or down to achieve decent acceptance rate
scale = 0.1 # higher for higher taus 
tuning = scale*basetuning
τ = 0.5
Σ = τ*(1.0-τ)*Z'Z/n
Σinv = inv(Σ)
θ = X\LNW  # OLS start values
# to do ordinary QR via MCMC, set Z=X (just to verify that MCMC works!)
lnL = θ -> likelihood(θ, LNW, X, Z, τ, Σinv)
Proposal = θ -> proposal(θ, tuning)
Prior = θ -> prior(θ)
chain = mcmc(θ, mcmcreps, burnin, Prior, lnL, Proposal)
d = dstats(chain)
# explore things:
p1 = plot(chain[:,1], legend=false, title="simple proposal: theta1") # chain is highly autocorrelated!
p2 = scatter(chain[:,1],chain[:,2], legend=false, title="simple proposal: theta1 vs. theta2") # there is significant posterior correlation
#= these results suggest that the proposal could be improved:
    * compute the covariance of the chain, and make draws from random walk MVN with that covariance
    * the following implements this idea
=#

# the following uses a different proposal, with correlations 
# across the parameters, based on the previous chain.
V = cov(chain[:,1:end-1])
cholV = cholesky(V).U
tuning = 0.5
Proposal = θ -> proposal2(θ, tuning*cholV)
chain = mcmc(θ, mcmcreps, burnin, Prior, lnL, Proposal)
p3 = plot(chain[:,1], legend=false, title="correlated proposal: theta1") # chain is highly autocorrelated!
p4 = scatter(chain[:,1],chain[:,2], legend=false, title="correlated proposal: theta1 vs. theta2") # there is significant posterior correlation
plot(p1,p2,p3,p4,layout=(2,2))
gui()
return
end
main()
