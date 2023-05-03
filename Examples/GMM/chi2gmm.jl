# GMM estimation for a sample from Chi^2(theta)
# compare to two method of moments estimators (see chi2mm.m)
using Distributions, Econometrics
n = 100
θ = 3.0
y = rand(Chisq(θ), n)
θhat = mean(y)
println("mm estimator based on mean: ", θhat)

θhat = 0.5*var(y)
println("mm estimator based on variance: ", θhat)

# GMM
W = eye(2)
obj = θ -> ([θ.-mean(y); θ.-0.5*var(y)]'W*[θ.-mean(y); θ.-0.5*var(y)]) 
θhat, junk, junk = fminunc(obj, [2.0])
println("GMM estimator based on mean and variance: ", θhat)

