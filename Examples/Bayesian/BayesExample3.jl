# sample model as in BayesExamples1 and 2, but using
# a specialized package and more sophisticated sampling
using Distributions, Turing, StatsPlots
function main()
# sample is from exponential, prior is lognormal
y = rand(Exponential(3.0), 30)
# the model: prior and likelihood
@model function ExpModel(y)
    θ ~ LogNormal(1.,1.)
    y ~ Exponential(θ)
end    
# get the chain
chain = sample(ExpModel(y), NUTS(200, 0.65), 10000)
end
## run it
chain = main()
# examine and plot chain
display(chain)
plot(chain[end-1000:end])

