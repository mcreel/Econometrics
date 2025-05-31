# sample model as in BayesExamples1 and 2, but using
# a specialized package and more sophisticated sampling

# note that we don't have to supply a proposal

##
using Distributions, Turing, StatsPlots
function main(n)
# sample is from exponential, prior is lognormal
y = rand(Exponential(3.0), n)
# the model: prior and likelihood
@model function ExpModel(y)
    θ ~ LogNormal(1.,1.) # the prior
    y ~ Exponential(θ)   # the likelihood
end    
# get the chain
chain = sample(ExpModel(y), NUTS(200, 0.65), 1000)
end

## run it
n = 300
chain = main(n)

## examine and plot last bit of chain chain
display(chain)
plot(chain[end-500:end])

