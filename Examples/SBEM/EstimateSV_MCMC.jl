using Econometrics, Plots, LinearAlgebra, Statistics, DelimitedFiles

include("SVlib.jl")

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

function main()
# setup
θtrue = TrueParameters() # true param values, on param space
lb, ub = PriorSupport()

# generate a new data set, or load the example data
n = 500
burnin = 100
#y = SVmodel(TrueParameters(),n, burnin) 
y = readdlm("svdata.txt")

# target statistics at sample data
m0 =  auxstat(y) 

# Estimation by indirect inference
S = 100 # number of simulation reps
shocks_u = randn(n+burnin,S) # used fixed shocks to control chatter
shocks_e = randn(n+burnin,S)
# define the moments
ms = θ -> m0' .- auxstat(θ, shocks_e, shocks_u)  # moment contributions: sample stats minus simulated stats
# do two-step GMM
# first stet
W = eye(size(m0,1))
m = θ -> vec(mean(ms(θ),dims=1)) # use average of moment contributions
obj = θ -> InSupport(θ) ? -m(θ)'*W*m(θ) : Inf

prior = θ -> Prior(θ)
# define things for MCMC
verbosity = true
burnin = 100
ChainLength = 1000
θinit = (ub + lb) ./ 2.0
# initial proposal moves one at a time
tuning = 0.2/sqrt(12.0)*(ub-lb) # two tenths of a standard. dev. of prior
proposal = θ -> proposal1(θ, tuning)
chain = mcmc(θinit, ChainLength, burnin, prior, obj, proposal, verbosity)
# now use a MVN random walk proposal 
Σ = NeweyWest(chain[:,1:3])
tuning = 1.0
for j = 1:5
    P = (cholesky(Σ)).U
    proposal = θ -> proposal2(θ,tuning*P)
    θinit = vec(mean(chain[:,1:3],dims=1))
    if j == 5
        ChainLength = 10000
    end    
    chain = mcmc(θinit, ChainLength, 0, prior, obj, proposal, verbosity)
    accept = mean(chain[:,end])
    #println("tuning: ", tuning, "  acceptance rate: ", accept)
    if accept > 0.35
        tuning *= 1.5
    elseif accept < 0.25
        tuning *= 0.5
    end
    θinit = vec(mean(chain[:,1:3],dims=1))
    Σ = 0.5*Σ + 0.5*cov(chain[:,1:3])
end
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
plot(p1,p2,p3)
plot!(legendfontsize=4)
#savefig("mcmc.png")
prettyprint([posmean inci])
return chain
end
main()
