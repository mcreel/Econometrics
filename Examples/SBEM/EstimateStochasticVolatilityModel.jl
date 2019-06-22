using Statistics, Calculus, LinearAlgebra, DelimitedFiles

include("SVlib.jl")

function main()
# setup
θtrue = [exp(-0.736/2.0), 0.9, 0.363] # true param values, on param space
lb = [0.0, 0.0, 0.0]
ub = [3.0, 1.0, 3.0]
n = 1000
burnin = 100
# load the example data
y = readdlm("svdata.txt")
m0 =  sqrt(n)*aux_stat(y) # statistic
# or use a new draw
#m0 = SVmodel(θtrue, n, randn(n+100), randn(n+100), true)   
# Estimation by indirect inference
S = 100 # number of simulation reps
shocks_u = randn(n+burnin,S)
shocks_e = randn(n+burnin,S)
ms = θ -> SVmoments(m0, n, θ, shocks_u, shocks_e)
weight = θ -> inv(cov(ms(θ)))
m = θ -> vec(mean(ms(θ),dims=1)) # 1Xg
obj = θ -> m(θ)'*weight(θ)*m(θ)
θhat, objvalue, converged, details = samin(obj, θtrue, lb, ub; ns = 5, verbosity = 2, rt = 0.25)
# compute the estimated standard errors and CIs
D = (Calculus.jacobian(m, vec(θhat), :central))
W = weight(θhat)
V = inv(D'*W*D)
se = sqrt.(diag(V))
println("true values, estimates, st. error, and limits of 95% CI")
prettyprint([θtrue θhat se θhat-1.96*se θhat+1.96*se],["true value", "estimate", "std. err.", "CI lower", "CI upper"])
end
main()
