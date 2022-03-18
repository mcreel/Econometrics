## This is an example of MSM that illustrates
#  the problem of "chatter" of the objective 
#  function, which makes extremum estimation
#  difficult. It goes on to explore Bayesian
#  MSM, which is pretty easy to implement.

#  We generate data that follows a Poisson
#  distribution, and estimate using some
#  simple simulated moments

# the standard parameterization of Poisson density
function λ(x,θ)
    exp.(x*θ)
end

## generate the "true" data
using Distributions
n = 100
x = randn(n,2)
θ₀ = [1.,1.] # true params
y = rand.(Poisson.(λ(x,θ₀)))

## the statistics, averaged over observations
function k(x,y)
    vec(mean([y sqrt.(y) x.*y x.*sqrt.(y)], dims=1))
end

## the simulated moments: S replications
function simulated_moments(θ, S, x)
    [k(x, rand.(Poisson.(λ(x,θ)))) for _ = 1:S]
end

## Define the MSM criterion 
S = 20 # number of simulation replications
m = θ -> k(x,y) - mean(simulated_moments(θ, S, x)) 
 # sums of squares of moments: corresponds to GMM with identity weight
function obj(θ)
    mean(abs2, m(θ))
end    
function obj(a,b) # a method for plots
    obj([a,b])
end

## profile the MSM objective function
using Plots
plotlyjs() # use this backend for interactive plots
θ₁ = range(0.9, length=100, stop=1.1)
θ₂ = θ₁
p1 = surface(θ₁, θ₂, (θ₁,θ₂)->obj(θ₁,θ₂),c=:viridis)
xlabel!("x")
ylabel!("y")
p2 = contour(θ₁, θ₂, (θ₁,θ₂)->obj(θ₁,θ₂),c=:viridis)
plot(p1, p2)

## The objective function looks horrible!
#  Nevertheless, let's attempt to do gradient-based estimation
using Econometrics
fminunc(obj, zeros(2))
# note that theta hat is same as start values, it didn't work

## Now, let's do Bayesian MSM, as suggested by
#  Chernozhukov and Hong (2003)
using Turing, AdvancedMH, Econometrics
k_data = k(x,y) # stats from the real data
#  Define the prior and asymptotic Gaussian likelihood of the moments
@model function MSM(k_data, S, x)
    θ ~ arraydist([Normal() for _=1:2]) # the prior (note: it's biased)
    # sample from the model, at the trial parameter value, and compute statistics
    ks = simulated_moments(θ, S, x)
    kbar = mean(ks) # mean of simulated statistics, over S replications
    Σ = cov(ks)  # estimated covariance of statistics
    k_data ~ MvNormal(kbar, Σ) # this corresponds to CUE GMM, but adding the determinant from the MVN 
end

# do MCMC sampling
# look at acceptance rate, below, and adjust
# the tuning to make acceptance around 0.25
S = 100 # make this as large as possible, given computation time, when doing real research
tuning = 0.01 
length = 10000
burnin = 1000
chain = sample(MSM(k_data, S, x),
    MH(:θ => AdvancedMH.RandomWalkProposal(MvNormal(zeros(2), tuning*eye(2)))),
    length; discard_initial=burnin)
display(chain)
plot(chain)
corner(chain)
c = Array(chain)
acceptance = size(unique(c[:,1]),1)[1] / size(c,1)
println("Acceptance rate: $acceptance")