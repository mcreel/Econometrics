# this computes the GMM estimator by SA minimization, for
# one of the 1000 data sets. Can be parallelized using
# threads.
using Econometrics, SolveDSGE, DelimitedFiles, Statistics, LinearAlgebra, ForwardDiff
include("DSGEmoments.jl")  # computes errors
include("CKlib.jl")

function main()

# using first line for the "true" data set, the second to generate new data
data = readdlm("dsgedata.txt")
#data = dgp(TrueParameters())

# define CUE GMM criterion
include("parameters.jl") # load true parameter values
lb, ub = PriorSupport()
n = size(data,1)
moments = theta -> DSGEmoments(theta, data)
m = theta -> vec(mean(moments(theta),dims=1)) # 1Xg
weight = theta -> inv(cov(moments(theta)))
obj = theta -> m(theta)'*weight(theta)*m(theta)
thetastart = (ub+lb)/2.0 # prior mean as start
# estimate by simulated annealing
thetahat, objvalue, converged, details = samin(obj, thetastart, lb, ub; ns = 20, verbosity = 1, rt = 0.5)
# compute the estimated standard errors and CIs
D = ForwardDiff.jacobian(m, vec(thetahat))
W = weight(thetahat)
V = inv(D'*W*D) 
se = sqrt.(diag(V))
println("estimates, st. error, and limits of 95% CI")
prettyprint([thetahat se thetahat-1.96*se thetahat+1.96*se],["estimate", "std. err.", "CI lower", "CI upper"])
# attempt gradient based (often doesn't work)
println("trying gradient-based")
thetahat, objvalue, converged = fmincon(obj, thetastart, [], [], lb, ub; iterlim=10000)
println("fmincon results. objective fn. value: ", objvalue)
println("parameter values:")
prettyprint(thetahat)
return nothing
end
main()
