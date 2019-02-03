""" 
    chain = mcmc(θ, reps, burnin, Prior, lnL, Proposal)

    Simple MH MCMC

    You must set the needed functions, e.g.,:
    Prior = θ -> your_prior(θ, whatever other args)
    lnL = θ -> your_log_likelihood(θ, whatever other args)
    Proposal = θ -> your_proposal(θ, whatever other args)
    (optionally) mcmcProposalDensity = (θᵗ,θ) -> your_proposal_density

    then get the chain using the above syntax, or optionally,
    (non-symmetric proposal) chain = mcmc(θ, reps, burnin, Prior, lnL, Proposal, ProposalDensity, report=true)
    (example code) mcmc(): runs a simple example. edit(mcmc,()) to see the code

"""
function mcmc()
    println("mcmc(), called with no arguments, runs a simple example")
    println("execute edit(mcmc,()) to see the code")
    # sample is from exponential, prior is lognormal, proposal is random walk lognormal
    y = rand(Exponential(3.0),30)
    # set prior, likelihood and proposal
    Prior = θ -> pdf.(Ref(LogNormal(1.0,1.0)), θ)
    lnL = θ -> sum(logpdf.(Ref(Exponential(θ)), y))
    tuning = 0.5
    Proposal = θ -> rand(LogNormal(log(θ),tuning))
    # get the chain, plot posterior, and descriptive stats
    chain = mcmc(1.0, 100000, 10000, Prior, lnL, Proposal, true) # start value, chain length, and burnin 
    p = npdensity(chain[:,1]) # nonparametric plot of posterior density 
    plot!(p, title="posterior density, simple MCMC example: true value = 3.0", show=true) # add a title
    dstats(chain)
    return
end

# method to allow skipping proposal density when symmetric
function mcmc(θ, reps, burnin, Prior, lnL, Proposal::Function, report=true::Bool)
    ProposalDensity = (a,b) -> 1.0
    mcmc(θ, reps, burnin, Prior, lnL, Proposal, ProposalDensity, report)
end    

# the main loop
function mcmc(θ, reps, burnin, Prior, lnL, Proposal::Function, ProposalDensity::Function, report=true::Bool)
    reportevery = Int((reps+burnin)/10)
    lnLθ = lnL(θ)
    chain = zeros(reps, size(θ,1)+1)
    naccept = zeros(size(θ))
    for rep = 1:reps+burnin
        θᵗ = Proposal(θ) # new trial value
        changed = Int.(.!(θᵗ .== θ)) # find which changed
        lnLθᵗ = lnL(θᵗ)
        # MH accept/reject
        accept = rand() < 
            exp(lnLθᵗ-lnLθ)*
            Prior(θᵗ)/Prior(θ)*
            ProposalDensity(θ,θᵗ)/ProposalDensity(θᵗ,θ)
        if accept
            θ = copy(θᵗ)
            lnLθ = lnLθᵗ 
        end
        naccept = naccept .+ changed .* Int.(accept)
        if (mod(rep,reportevery)==0 && report)
            println("current parameters: ", round.(θ,digits=3))
            println("  acceptance rates: ", round.(naccept/reportevery,digits=3))
            naccept = naccept - naccept
        end    
        if rep > burnin
            chain[rep-burnin,:] = [θ; accept]
        end    
    end
    return chain
end


