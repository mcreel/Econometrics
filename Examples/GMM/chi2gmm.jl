# GMM estimation for a sample from Chi^2(theta)
# compare to two method of moments estimators (see chi2mm.m)
using Distributions
n = 30
theta = 3.0
y = rand(Chisq(theta), n)
thetahat = mean(y)
println("mm estimator based on mean: ", thetahat)

thetahat = 0.5*var(y)
println("mm estimator based on variance: ", thetahat)

# GMM
W = eye(2)
obj = theta -> ([theta.-mean(y); theta.-0.5*var(y)]'W*[theta.-mean(y); theta.-0.5*var(y)]) 
thetahat, junk, junk = fminunc(obj, [2.0])
println("GMM estimator based on mean and variance: ", thetahat)

