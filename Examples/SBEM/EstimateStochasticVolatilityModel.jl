using Statistics, Calculus, LinearAlgebra, DelimitedFiles

include("SVlib.jl")

#function main()
# generate the sample
θtrue = [exp(-0.736/2.0), 0.9, 0.363] # true param values, on param space
lb = [0.0, 0.0, 0.0]
ub = [3.0, 1.0, 3.0]
burnin = 100
# load the example data
y = readdlm("svdata.txt")
y = y[:]
n = size(y,1) # sample size
m0 =  aux_stat(y) # generate the sample and save the data
# Estimation by indirect inference
S = 100 # number of simulation reps
randdraws = randn(S,n+burnin,2) # fix the shocks to control "chatter" (includes the burnin period)
ms = θ -> SVmoments(m0, n, θ, randdraws)
m = θ -> vec(mean(ms(θ),dims=1)) # 1Xg
weight = θ -> inv(cov(ms(θ)))
obj = θ -> m(θ)'weight(θ)*m(θ)
thetahat, objvalue, converged, details = samin(obj, θtrue, lb, ub; ns = 5, verbosity = 2, rt = 0.5)
# compute the estimated standard errors and CIs
D = (Calculus.jacobian(m, vec(thetahat), :central))
W = weight(thetahat)
V = inv(D'*W*D)
se = sqrt.(diag(V))

println("true values, estimates, st. error, and limits of 95% CI")
prettyprint([θtrue thetahat se thetahat-1.96*se thetahat+1.96*se],["true value", "estimate", "std. err.", "CI lower", "CI upper"])
#end
#main()
