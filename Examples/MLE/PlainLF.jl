# this function shows the behavior of the 
# likelihood, the log-likelihood and the average
# log-likelihood as the sample size increases
using Distributions
n = 10 # try increasing this to 50, 1000, etc.
p = 0.3
y = rand(Bernoulli(p),n)
f =  pdf.(Bernoulli(p), y)
@show n
@show ℒ = prod(f)
@show lnℒ = sum(log.(f)) 
@show avglnℒ = (1/n)*sum(log.(f))
nothing


