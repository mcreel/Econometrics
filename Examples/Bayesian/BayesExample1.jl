# Example of Bayesian estimation
# sampling from exponential(theta) distribution
# lognormal prior

# Shows how likelihood, joint, marginal likelihood,
# and posterior are formed, and how posterior can be
# integrated to get posterior mean

# explore different sample sizes, different true thetas
using Distributions, Plots

function main()
    n = 50   # sample size
    theta = 3 # true theta
    y = rand(Exponential(theta), n) # sample from exponential(theta)
    # make plots
    thetas = range(0.01,stop=10,length=1000)
    delta = thetas[2]-thetas[1]
    p, priormean = prior(thetas)
    post, postmean = posterior(y, thetas)
    plot(thetas, [p post], label = ["prior" "posterior"])
    plot!([priormean, priormean], [0.0, 1.0], label = "prior mean")
    plot!([postmean, postmean], [0.0, 1.0], label = "posterior mean")

end

# the prior is lognormal(1,1)
function prior(theta)
    d = LogNormal(1.0,1.0)
    p = pdf.(Ref(d), theta)
    pmean = exp(1.5) # mean of lognormal is exp(mu+sig/2)
    return p, pmean
end

# the likelihood function
function likelihood(y, theta)
    dens = zeros(size(theta))
    for i = 1:size(theta,1)
        d = Exponential(theta[i])
        dens[i] = prod(pdf.(Ref(d), y))
    end
    return dens
end


# joint is prior X likelihood
function joint(y, theta)
    l = likelihood(y, theta)
    p, junk = prior(theta)
    dens = l.*p
end

# compute marginal likelihood of Y by integrating out theta (crude, only illustrative)
function marginal(y)
    dens = 0.0
    theta = 0.0
    delta = 0.01
    # evaluate joint over grid
    for r = 1:1000
        theta += delta
        dens += joint(y, theta)
    end
    # marginalize by integrating the joint (sum up height X width)
    dens = dens*delta
end

# the posterior, by Bayes' Law
function posterior(y, theta)
    m = marginal(y)
    j = joint(y, theta)
    dens = j ./ m
    thetas = range(0.01,stop=10,length=1000)
    pmean = sum(dens.*thetas.*0.01)
    return dens, pmean
end

main()

