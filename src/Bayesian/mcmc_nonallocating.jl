# this is an experimental version that allocates less. It is 
# about 30% faster, for a very simple likelihood, but the
# time is so short that it's not really important. For more complex
# likelihoods, it seems likely that making the likelihood computation fast is more beneficial than making the sampler fast (?)
using Distributions

function Prior!(p, θ)
    p .= pdf.(LogNormal(1.0,1.0), θ[1])
end

function lnL!(ll, θ, y)
    ll .= sum(logpdf.(Exponential(θ[1]), y))
end

function Proposal!(θᵗ, θᶜ, tuning)
    θᵗ .= rand(LogNormal(log(θᶜ[1]),tuning))
end

struct model
    data::Array{Float64}
    prior!::Function
    lnL!::Function
    proposal!::Function
end    

function mcmc(θᶜ, model, reps, burnin)
    tuning = 0.5
    θᵗ = copy(θᶜ) # θᵗ params
    chain = copy(θᶜ) # need defined at this scope
    ℒᵗ = [-Inf]
    ℒᶜ = [-Inf]
    pt = [-Inf]
    pc = [1.]
    Prior!(pc, θᶜ)
    model.lnL!(ℒᶜ, θᶜ, model.data) # get initial likelihood
    for i = 1:reps + burnin
        model.proposal!(θᵗ, θᶜ, tuning)
        model.lnL!(ℒᵗ, θᵗ, model.data)
        Prior!(pt, θᵗ)
        # if accepted, update the current quantities
        if rand() < exp(ℒᵗ[1] - ℒᶜ[1])*pt[1]/pc[1]
            θᶜ .= θᵗ
            pc .= pt 
            ℒᶜ .= ℒᵗ
        end
        # if out of burning, push current to chain
        i > burnin + 1 ? push!(chain, θᶜ[1]) : 
            (i == burnin + 1 ? chain = copy(θᶜ) : nothing)
    end
    acceptance = (size(unique(chain)) ./ size(chain))[1]
return chain, acceptance
end

function main()
y = rand(Exponential(3.0),30)
θ = [rand(LogNormal(1.0, 1.0))]   # initial value draw from prior
m = model(y, Prior!, lnL!, Proposal!)
chain, acceptance = mcmc(θ, m, 100000, 10000)
end

