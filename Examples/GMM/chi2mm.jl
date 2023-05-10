# method of moments estimators for a sample from Chi^2(theta)
# increase n to observe consistency
# note that the two estimators are different from one another
##
using Distributions, StatsPlots
n = 100
theta = 3

# the distribution doesn't matter, what matters is
# that the moments are correctly specified
#y = rand(Normal(theta, sqrt(2*theta)), n) # this Normal function takes μ and σ
y = rand(Chisq(theta), n)
display(density(y))

# a MM estimator
thetahat = mean(y)
println("mm estimator based on mean: ", thetahat)
vline!([thetahat])

# a second MM estimator
thetahat = 0.5*var(y)
println("mm estimator based on variance: ", thetahat)
vline!([thetahat])
