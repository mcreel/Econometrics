# requires the package https://github.com/mcreel/SV
using SV, Statistics, Calculus, LinearAlgebra, DelimitedFiles

function main()
# setup
θtrue = [exp(-0.736/2.0), 0.9, 0.363] # true param values, on param space
lb = [0.0, 0.0, 0.0]
ub = [2.0, 0.99, 1.0]
n = 1000
burnin = 100
# load the example data
y = readdlm("svdata.txt")
m0 =  sqrt(n)*aux_stat(y) # statistic
@show m0
# or use a new draw
#m0 = SVmodel(θtrue, n, randn(n+100), randn(n+100), true)   
# Estimation by indirect inference
S = 100 # number of simulation reps
shocks_u = randn(n+burnin,S)
shocks_e = randn(n+burnin,S)
ms = θ -> SVmoments(m0, n, θ, shocks_u, shocks_e)
# the estimation of the covariance is unstable when volatility is 
# high, so first, identity weight, then efficient, with tighter bounds
weight = θ -> 1.0
m = θ -> vec(mean(ms(θ),dims=1)) # 1Xg
obj = θ -> m(θ)'*weight(θ)*m(θ)
θinit = (lb + ub)./2.0
θhat, objvalue, converged, details = samin(obj, θinit, lb, ub; ns = 5, verbosity = 3, rt = 0.5)
lb2 = 0.9.*θhat
ub2 = 1.1.*θhat
lb = max.(lb,lb2)
ub = min.(ub, ub2)
weight = θ -> inv(cov(ms(θ)))
θhat, objvalue, converged, details = samin(obj, θhat, lb, ub; ns = 5, verbosity = 3, rt = 0.5)
# compute the estimated standard errors and CIs
D = (Calculus.jacobian(m, vec(θhat), :central))
W = weight(θhat)
V = inv(D'*W*D)
se = sqrt.(diag(V))
println("true values, estimates, st. error, and limits of 95% CI")
prettyprint([θtrue θhat se θhat-1.96*se θhat+1.96*se],["true value", "estimate", "std. err.", "CI lower", "CI upper"])
end
main()
