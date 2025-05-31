## this computes the GMM estimator by SA minimization, using
# the default typical data set (or generating new data)
using Econometrics, SolveDSGE, CSV, Statistics, LinearAlgebra, ForwardDiff
cd(@__DIR__)
include("DSGEmoments.jl")  # computes errors
include("CKlib.jl")
# needed to explore other data sets 
global const dsge = retrieve_processed_model("CK_processed.txt")

function main(defaultdata = true)
cd(@__DIR__)
# call program as main() for the default data, as main(false) for new random data
defaultdata ? data = CSV.File("dsgedata.csv") |> CSV.Tables.matrix : data = dgp(TrueParameters(), dsge, 1, rand(1:Int64(1e10)))[1]

# define CUE GMM criterion
include("parameters.jl") # load true parameter values
lb, ub = PriorSupport()
n = size(data,1)
moments = theta -> sqrt(n)*DSGEmoments(theta, data)
m = theta -> vec(mean(moments(theta),dims=1)) # 1Xg
weight = theta -> inv(cov(moments(theta)))
obj = theta -> m(theta)'*weight(theta)*m(theta)
thetastart = (ub+lb)/2.0 # prior mean as start

# estimate by simulated annealing
thetahat, objvalue, converged, details = samin(obj, thetastart, lb, ub; ns = 20, nt=5, verbosity = 1, rt = 0.75)
# compute the estimated standard errors and CIs
D = ForwardDiff.jacobian(m, vec(thetahat))
W = weight(thetahat)
V = inv(D'*W*D)/n 
se = sqrt.(diag(V))
println("estimates, st. error, and limits of 95% CI")
prettyprint([thetahat se thetahat-1.96*se thetahat+1.96*se],["estimate", "std. err.", "CI lower", "CI upper"])
# attempt gradient based (often doesn't work)
println("trying gradient-based")
thetahat2, objvalue, converged = fmincon(obj, thetastart, [], [], lb, ub; iterlim=10000)
println("fmincon results. objective fn. value: ", objvalue)
println("parameter values:")
prettyprint(thetahat2)
return thetahat
end
θhat = main();
