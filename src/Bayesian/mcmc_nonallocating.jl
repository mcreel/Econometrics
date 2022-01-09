# this is an experimental version that allocates less. It is 
# about twice as fast, for a very simple likelihood, but the
# time is so short that it's not really beneficial. It seems
# likely that making the likelihood computation fast is more
# important.
using Distributions

function Prior!(p, θ)
    p .= pdf.(LogNormal(1.0,1.0), θ[1])
end

function lnL!(ll, θ, y)
    ll .= sum(logpdf.(Exponential(θ[1]), y))
end

function Proposal!(trial, current, tuning)
    trial .= rand(LogNormal(log(current[1]),tuning))
end

struct model
    data::Array{Float64}
    prior!::Function
    lnL!::Function
    proposal!::Function
end    

function mcmc(θ, model, reps, burnin)
    tuning = 0.5
    current = copy(θ)
    trial = copy(θ)
    chain = copy(θ)
    lt = [-Inf]
    lc = [-Inf]
    model.lnL!(lc, current, model.data) # get initial likelihood
    for i = 1:reps + burnin
        model.proposal!(trial, current, tuning)
        model.lnL!(lt, trial, model.data)
        if rand() < exp(lt[1] - lc[1])   #### NEED TO ADD PRIOR RATIO
            current .= trial
            lc .= lt
        end
        i == burnin + 1 ? chain = copy(current) : nothing
        i > burnin + 1 ? push!(chain, current[1]) : nothing   
    end
return chain
end

function main()
y = rand(Exponential(3.0),30)
θ = [rand(LogNormal(1.0, 1.0))]   # initial value draw from prior
m = model(y, Prior!, lnL!, Proposal!)
chain = mcmc(θ, m, 100000, 10000)
#display(chain)
end
# the main loop

