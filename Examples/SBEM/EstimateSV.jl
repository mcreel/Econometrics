using Econometrics, Statistics, Calculus, LinearAlgebra, DelimitedFiles

include("SVlib.jl") # library of functions for the SV model

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
obj = θ -> InSupport(θ) ? m(θ)'*W*m(θ) : Inf
θinit = (lb + ub)./2.0 
θhat, objvalue, converged, details = samin(obj, θinit, lb, ub; ns = 5, verbosity = 3, rt = 0.5)
# second step 
W = inv(cov(ms(θhat))) # get the optimal weight matrix
obj = θ -> InSupport(θ) ? m(θ)'*W*m(θ) : Inf
θhat, objvalue, converged, details = samin(obj, θhat, lb, ub; ns = 5, verbosity = 3, rt = 0.5)
# compute the estimated standard errors and CIs
D = (Calculus.jacobian(m, vec(θhat), :central))
V = inv(D'*W*D)
se = sqrt.(diag(V))
println("true values, estimates, st. error, and limits of 95% CI")
prettyprint([θtrue θhat se θhat-1.96*se θhat+1.96*se],["true value", "estimate", "std. err.", "CI lower", "CI upper"])
end
main()
