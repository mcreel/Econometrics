# this function shows the behavior of the 
# likelihood, the log-likelihood and the average
# log-likelihood as the sample size increases
using Distributions
function main(n)
p = 0.3
y = rand(Bernoulli(p),n)
f =  pdf(Bernoulli(p), y)
@show lnℒ = sum(log.(f)) 
@show avglnℒ = (1/n)*sum(log.(f)) 
@show ℒ = prod(f)
nothing
end


