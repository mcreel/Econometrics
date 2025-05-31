## This provides examples for MSM and SML.
#  The first is an example for MSM that illustrates
#  the problem of "chatter" of the objective 
#  function, which makes extremum estimation
#  difficult. It goes on to explore Bayesian
#  MSM, which is pretty easy to implement.
#  
#  The second is for SML, and shows how simulated
#  discrete dependent variable models can suffer 
#  from discontinuities in the criterion. This can also
#  be dealt with using the Bayesian approach.


########################### MSM Example ###############################
using Econometrics, Term, Pkg
cd(@__DIR__)
Pkg.activate("../")
println(@yellow "Start of MSM example. We generate data that follows a Gaussian distribution, and estimate using some simple simulated moments")


## DATA
using Distributions
println(@yellow "Generate the 'true' data. The conditional mean is usual linear function")
function μ(x,θ)
    x*θ
end
n = 100
x = randn(n,2)
θ₀ = [1.,1.] # true params
y = rand.(Normal.(μ(x,θ₀),1.0))

## STATISTICS 
# the statistics, averaged over observations: these are just off
# the top of the head, they may not be very good
function Z(x,y)
    vec(mean([y x.*y x.*(y.^2.0) (x.^2.0).*y], dims=1))
end
# the simulated moments: S replications
using Random
function simulated_moments(θ, x, S, controlchatter)
    controlchatter ? Random.seed!(1234) : nothing
    [Z(x, rand.(Normal.(μ(x,θ),1.0))) for _ = 1:S]
end

## MSM SETUP 
# Define the and profile MSM criterion 
# IMPORTANT: try running this with chatter, and without.
S = 100 # number of simulation replications
controlchatter = false # try setting this to false/true, and run this and the next block
m = θ -> Z(x,y) - mean(simulated_moments(θ, x, S, controlchatter)) 
 # sums of squares of moments: corresponds to GMM with identity weight
function obj(θ)
    mean(abs2, m(θ))
end    
function obj(a,b) # a method for plots
    obj([a,b])
end
# profile the MSM objective function
using Plots
gr()
θ₁ = range(0., length=100, stop=2.0)
θ₂ = θ₁
p2 = contour(θ₁, θ₂, (θ₁,θ₂)->obj(θ₁,θ₂),c=:viridis)
xlabel!("x")
ylabel!("y")


## ESTIMATION
# The objective function looks horrible, when there is chatter
# without chatter, it looks good. Experiment estimation with and without
# chatter, bu
#  Nevertheless, let's attempt to do gradient-based estimation
using Econometrics
# run this next line repeatedly, using CTRL-enter. When chatter
# is controlled, you always get the same result. When it's not,
# the result changes every time.
fminunc(obj, ones(2))
println(@green "You should find that MSM works, when chatter is controlled")
println(@yellow "(Remember that the true parameters are $θ₀)")
# note that, with chatter, the estimates change every time, even though the sample does not.

## Now, let's do Bayesian MSM, as suggested by
#  Chernozhukov and Hong (2003)
# This works even without controlling chatter
using Turing, AdvancedMH, Term, Econometrics
Z_data = Z(x,y) # stats from the real data
#  Define the prior and approximate Gaussian likelihood of the moments
@model function MSM(Z_data, S, x)
    θ ~ arraydist([Normal() for _=1:2]) # the prior (note: it's biased)
    # sample from the model, at the trial parameter value, and compute statistics
    Zs = simulated_moments(θ, x, S, false) # NOTE: no control of chatter
    Zbar = mean(Zs) # mean of simulated statistics, over S replications
    Σ = cov(Zs)  # estimated covariance of statistic
    Z_data ~ MvNormal(Zbar, Σ) # this corresponds to CUE GMM, but adding the determinant from the MVN 
end
# do MCMC sampling
# look at acceptance rate, below, and adjust
# the tuning to make acceptance around 0.25
S = 50 # make this as large as possible, given computation time, when doing real research
tuning = 0.05
length = 2000
burnin = 200
chain = sample(MSM(Z_data, S, x),
    MH(:θ => AdvancedMH.RandomWalkProposal(MvNormal(zeros(2), tuning*eye(2)))),
    length; init_params=zeros(2), discard_initial=burnin)
display(chain)
plot(chain)
#corner(chain)
c = Array(chain)
acceptance = size(unique(c[:,1]),1)[1] / size(c,1)
println("Acceptance rate: $acceptance")
println(@red "you want to see acceptance around 0.2, and ESS above 200")
println(@red "Increase S and adjust tuning to achieve that. Try S=100, tuning=0.05")
println(@cyan "End of MSM example")

########################### SML Example ###############################

## Now an SML example, using a simple Logit model,
#  where there is no actual need for simulation. 
#  However, it does show how simulation can introduce
#  discontinuity. This example controls random numbers,
#  by re-setting the seed for each simulation, which 
#  eliminates the chatter problem seen above. However,
#  discontinuity still causes problems for extremum 
#  estimators.

## generate the "true" data
using Term
println(@yellow "Start of SML example")
using Distributions, Random
n = 200
x = randn(n,2)
θ₀ = [1.,1.] # true params
# get the responses
response(θ,x) = rand(size(x,1)) .< 1.0./(1.0 .+ exp.(-x*θ)) 
y = response(θ₀,x) 

## define functions for SML
# define simulator of probabilities
simulated_probs(θ, x, S) = mean([response(θ, x) for _ = 1:S])
# define the simulated log likelihood
# IMPORTANT: note how the seed is re-set before simulating
# the choice probabilities. This eliminates chatter.
function logit_simulated_loglikelihood(θ, y, x, S)
    Random.seed!(1234)
    p = simulated_probs(θ, x, S) 
    y.*log.(p) .+ (log.(1.0 .- p)).*(1.0 .- y)
end
# define the objective function for SML (minimization form)
S = 100
function smlobj(θ)
    -mean(logit_simulated_loglikelihood(θ, y, x, S))
end    
function smlobj(a,b) # method for plots
    smlobj([a,b])
end    

## profile the SML objective function
using Plots
gr()
θ₁ = range(0.6, length=100, stop=1.2)
θ₂ = θ₁
p2 = contour(θ₁, θ₂, (θ₁,θ₂)->smlobj(θ₁,θ₂),c=:viridis)
xlabel!("x")
ylabel!("y")
## The objective function looks horrible!
#  Nevertheless, let's attempt to do gradient-based estimation
using Econometrics
display(fminunc(smlobj, zeros(2)))
println(@red "Bzzzt! That didn't work!")

## Now, let's do Bayesian SML, again using 
#  methods from Chernozhukov and Hong (2003)
using Turing, AdvancedMH
#  Define the prior and asymptotic Gaussian likelihood of the moments
@model function SML(y, x, S)
    θ ~ arraydist([Normal() for _=1:2]) # the prior (note: it's biased)
    p = simulated_probs(θ, x, S)
    for i in eachindex(y)
        y[i] ~ Bernoulli(p[i])
    end
end
# do MCMC sampling
# look at acceptance rate, below, and adjust
# the tuning to make acceptance around 0.25
S = 200 # make this as large as possible, given computation time, when doing real research
tuning = 0.1
length = 2500
burnin = 500
chain = sample(SML(y, x, S),
    MH(:θ => AdvancedMH.RandomWalkProposal(MvNormal(zeros(2), tuning*eye(2)))),
    length; discard_initial=burnin)
display(chain)
plot(chain)
#corner(chain)
c = Array(chain)
acceptance = size(unique(c[:,1]),1)[1] / size(c,1)
println("Acceptance rate: $acceptance")
println(@cyan "End of SML example")
