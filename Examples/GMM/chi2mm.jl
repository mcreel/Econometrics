# method of moments estimators for a sample from Chi^2(theta)
# increase n to observe consistency
# note that the two estimators are different from one another
using Distributions
n = 30
theta = 3

# the distribution doesn't matter, what matters is
# that the moments are correctly specified
y = rand(Chisq(theta), n)
#y = rand(Normal(theta, sqrt(2*theta)), n)

# a MM estimator
thetahat = mean(y)
println("mm estimator based on mean: ", thetahat)

# a second MM estimator
thetahat = 0.5*var(y)
println("mm estimator based on variance: ", thetahat)
