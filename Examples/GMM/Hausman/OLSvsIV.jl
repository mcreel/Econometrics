# this is a little Monte Carlo exercise that illustrates that
# the OLS estimator is biased and inconsistent when errors are
# correlated with regressors, but that the IV estimator is consistent
using Econometrics, Plots, LinearAlgebra
function main()
reps = 1000 # number of Monte Carlo reps.
betaols = zeros(reps,2)
betaiv = zeros(reps,2)
n = 1000 # sample size

# covariance of X, W, e
cor_X_W = 0.2  # experiment with lowering or raising this: quality of instrument
cor_X_e = 0.5  # try setting this to zero to observe inefficiency of IV
sig = [
      1.0       cor_X_W     cor_X_e;
      cor_X_W   1.0         0.0;
     cor_X_e    0.0         1.0]
truebeta = [1, 2] # true beta
p = cholesky(sig).U
for i = 1:reps
	XWE = randn(n,3)*p
	e = XWE[:,3:3]
	w = [ones(n) XWE[:,2:2]]
	x = [ones(n) XWE[:,1:1]]
	y = x*truebeta + e
	# OLS
	betaols[i,:] = (x\y)'
	# IV
	betaiv[i,:] = (inv(w'*x)*w'*y)'
end

p1 = npdensity(betaols[:,2])
plot!(p1,title="OLS")
p2 = npdensity(betaiv[:,2])
plot!(p2,title="IV")
p = plot(p1,p2,layout=(2,1))
#savefig("olsiv.png")
gui()
println("true betas are ", truebeta)
println("OLS results")
dstats(betaols, short=true)
println("IV results")
dstats(betaiv, short=true)
return p
end
main()
