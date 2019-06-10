# this computes the GMM estimator by SA minimization, for
# one of the 1000 data sets. Can be parallelized using
# threads.
using Econometrics, DelimitedFiles, Statistics, LinearAlgebra, Calculus
include("DSGEmoments.jl")  # computes errors
function main()
include("parameters.jl") # load true parameter values
data = readdlm("dsgedata.txt")
n = size(data,1)
# estimate by simulated annealing
lb = lb_param_ub[:,1]
ub = lb_param_ub[:,3]
# define CUE GMM criterion
moments = theta -> DSGEmoments(theta, data)
m = theta -> vec(mean(moments(theta),dims=1)) # 1Xg
weight = theta -> inv(cov(sqrt(n)*moments(theta)))
obj = theta -> m(theta)'*weight(theta)*m(theta)
thetastart = (ub+lb)/2.0 # prior mean as start
# attempt gradient based (doesn't work)
thetahat, objvalue, converged = fmincon(obj, thetastart, [], [], lb, ub; iterlim=10000)
println("fmincon results. objective fn. value: ", objvalue)
println("parameter values:")
prettyprint(thetahat)
# simulated annealing
thetahat, objvalue, converged, details = samin(obj, thetastart, lb, ub; ns = 20, verbosity = 2, rt = 0.5)
# compute the estimated standard errors and CIs
D = (Calculus.jacobian(m, vec(thetahat), :central))
W = weight(thetahat)
V = inv(D'*W*D)/n 
se = sqrt.(diag(V))
println("estimates, st. error, and limits of 95% CI")
prettyprint([thetahat se thetahat-1.96*se thetahat+1.96*se],["estimate", "std. err.", "CI lower", "CI upper"])
return nothing
end
main()
