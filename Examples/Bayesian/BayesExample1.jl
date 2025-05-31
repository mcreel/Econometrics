# Example of Bayesian estimation
# sampling from exponential(θ) distribution
# lognormal prior

# Shows how likelihood, joint, marginal likelihood,
# and posterior are formed, and how posterior can be
# integrated to get posterior mean

## explore different sample sizes, different true θs
using Econometrics, Distributions, Plots

## the prior is lognormal(1,1)
function prior(θ)
    μ = 1.0
    σ = 1.0
    d = LogNormal(μ,σ)
    p = pdf.(Ref(d), θ)
    pmean = exp(μ + σ/2.0) # mean of lognormal is exp(mu+sig/2)
    return p, pmean
end

# the likelihood function
function likelihood(y, θ)
    dens = zeros(size(θ))
    for i ∈ axes(θ,1)
        d = Exponential(θ[i])
        dens[i] = prod(pdf.(Ref(d), y))
    end
    return dens
end

# joint is prior X likelihood
function joint(y, θ)
    l = likelihood(y, θ)
    p, junk = prior(θ)
    dens = l.*p
end

# compute marginal likelihood of Y by integrating out θ (crude, only illustrative)
function marginal(y)
    dens = 0.0
    θ = 0.0
    delta = 0.01
    # evaluate joint over grid
    for r = 1:1000
        θ += delta
        dens += joint(y, θ)  # doing the summing here (see 3 lines down)
    end
    # marginalize by integrating the joint (sum up height X width)
    dens = dens*delta        # rather than here
end

# the posterior, by Bayes' Law
function posterior(y, θ)
    m = marginal(y)
    j = joint(y, θ)
    dens = j ./ m       # get the conditional density
    θs = range(0.01,stop=10,length=1000)
    pmean = sum(dens.*θs.*0.01)  # crude sum rectangle integration
    return dens, pmean
end

function main()
    n = 100   # sample size
    θ = 3.0 # true θ
    y = rand(Exponential(θ), n) # sample from exponential(θ)
    # make plots
    θs = range(0.01,stop=10,length=1000)
    p, priormean = prior(θs)
    post, postmean = posterior(y, θs)
    plot(θs, [p post], label = ["prior" "posterior"])
    plot!([priormean, priormean], [0.0, 1.0], label = "prior mean")
    plot!([postmean, postmean], [0.0, 1.0], label = "posterior mean")
    plot!([θ, θ], [0.0, 1.0], label = "θ₀")

end

main()

