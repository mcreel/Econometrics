using Econometrics, Plots, Statistics, Calculus, LinearAlgebra, DelimitedFiles

include("SVlib.jl") # library of functions for the SV model

# proposals for MCMC
# uniform random walk in one dimension
function proposal1(current, tuning)
    trial = copy(current)
    i = rand(1:size(trial,1))
    trial[i] += tuning[i].*randn()
    return trial
end
# MVN random walk, or occasional draw from prior
function proposal2(current, cholV)
    current + cholV'*randn(size(current))
end

function main()
# setup
θtrue = TrueParameters() # true param values, on param space
lb, ub = PriorSupport()

# generate a new data set, or load the example data
n = 500
burnin = 100
mcreps = 100
results = zeros(mcreps,22)
Threads.@threads for rep = 1:mcreps
    y = SVmodel(TrueParameters(),n, burnin) 
    m0 =  auxstat(y) 
    # Estimation by two-step MSM
    S = 10 # number of simulation reps
    shocks_u = randn(n+burnin,S) # used fixed shocks to control chatter
    shocks_e = randn(n+burnin,S)
    # define the moments
    ms = θ -> m0' .- auxstat(θ, shocks_e, shocks_u)  # moment contributions: sample stats minus simulated stats
    # first step
    W = eye(size(m0,1))
    m = θ -> vec(mean(ms(θ),dims=1)) # use average of moment contributions
    obj = θ -> InSupport(θ) ? m(θ)'*W*m(θ) : Inf
    θhat, objvalue, converged, details = samin(obj, θtrue, lb, ub; ns = 5, verbosity = 0, rt = 0.5)
    # Estimate the weight matrix using more draws
    W = inv((1.0 + 1.0/S)*cov(m0' .- auxstat(θhat, randn(n+burnin, 500), randn(n+burnin, 500))))
    # second step
    obj = θ -> InSupport(θ) ? m(θ)'*W*m(θ) : Inf
    θhat, objvalue, converged, details = samin(obj, θhat, lb, ub; ns = 5, verbosity = 0, rt = 0.5)
    # compute the estimated standard errors and CIs
    D = (Calculus.jacobian(m, vec(θhat), :central))
    se = sqrt.(diag(inv(D'*W*D)))
    # now do MCMC
    prior = θ -> Prior(θ)
    burnin = 100
    ChainLength = 1000
    # initial proposal moves one at a time
    tuning = 0.2/sqrt(12.0)*(ub-lb) # two tenths of a standard. dev. of prior
    proposal = θ -> proposal1(θ, tuning)
    obj = θ -> InSupport(θ) ? -m(θ)'*W*m(θ) : -Inf
    chain = mcmc(θhat, ChainLength, burnin, prior, obj, proposal, false)
    # now use a MVN random walk proposal 
    Σ = NeweyWest(chain[:,1:3])
    tuning = 1.0
    mcmcreps = 3
    for j = 1:mcmcreps
        P = (cholesky(Σ)).U
        proposal = θ -> proposal2(θ,tuning*P)
        θinit = vec(mean(chain[:,1:3],dims=1))
        if j == mcreps
            ChainLength = 10000
        end    
        chain = mcmc(θinit, ChainLength, 0, prior, obj, proposal, false)
        if j < mcreps
            accept = mean(chain[:,end])
            if accept > 0.35
                tuning *= 1.5
            elseif accept < 0.25
                tuning *= 0.5
            end
            θinit = vec(mean(chain[:,1:3],dims=1))
            Σ = 0.5*Σ + 0.5*cov(chain[:,1:3])
        end    
    end
    q1 = quantile(chain[:,1],[0.005, 0.025, 0.05, 0.95, 0.975, 0.995])
    q2 = quantile(chain[:,2],[0.005, 0.025, 0.05, 0.95, 0.975, 0.995])
    q3 = quantile(chain[:,3],[0.005, 0.025, 0.05, 0.95, 0.975, 0.995])
    bayes_inci = [
            θtrue[1] >= q1[1] && θtrue[1] <= q1[6],
            θtrue[2] >= q2[1] && θtrue[2] <= q2[6],
            θtrue[3] >= q3[1] && θtrue[3] <= q3[6],
            θtrue[1] >= q1[2] && θtrue[1] <= q1[5],
            θtrue[2] >= q2[2] && θtrue[2] <= q2[5],
            θtrue[3] >= q3[2] && θtrue[3] <= q3[5],
            θtrue[1] >= q1[3] && θtrue[1] <= q1[4],
            θtrue[2] >= q2[3] && θtrue[2] <= q2[4],
            θtrue[3] >= q3[3] && θtrue[3] <= q3[4]
           ]
    extremum_inci = [
            θtrue[1] >= θhat[1] - 2.576*se[1] && θtrue[1] <= θhat[1] + 2.576*se[1],                
            θtrue[2] >= θhat[2] - 2.576*se[2] && θtrue[2] <= θhat[2] + 2.576*se[2],                
            θtrue[3] >= θhat[3] - 2.576*se[3] && θtrue[3] <= θhat[3] + 2.576*se[3],                
            θtrue[1] >= θhat[1] - 1.96 *se[1] && θtrue[1] <= θhat[1] + 1.96* se[1],                
            θtrue[2] >= θhat[2] - 1.96* se[2] && θtrue[2] <= θhat[2] + 1.96* se[2],                
            θtrue[3] >= θhat[3] - 1.96* se[3] && θtrue[3] <= θhat[3] + 1.96* se[3],                
            θtrue[1] >= θhat[1] - 1.645*se[1] && θtrue[1] <= θhat[1] + 1.645*se[1],                
            θtrue[2] >= θhat[2] - 1.645*se[2] && θtrue[2] <= θhat[2] + 1.645*se[2],                
            θtrue[3] >= θhat[3] - 1.645*se[3] && θtrue[3] <= θhat[3] + 1.645*se[3]
            ]
   results[rep,:] = vcat(converged, θhat, extremum_inci, bayes_inci)
end
return results
end
results = main()
extremum_inci = mean(results[:,5:13],dims=1)
bayes_inci = mean(results[:,14:22],dims=1)
println("99% CI, extremum")
prettyprint(extremum_inci[1:3], ["ϕ"; "ρ"; "σ"])
println("99% CI, Bayes")
prettyprint(bayes_inci[1:3], ["ϕ"; "ρ"; "σ"])
println("95% CI, extremum")
prettyprint(extremum_inci[4:6], ["ϕ"; "ρ"; "σ"]) 
println("95% CI, Bayes")
prettyprint(bayes_inci[4:6], ["ϕ"; "ρ"; "σ"])
println("90% CI, extremum")
prettyprint(extremum_inci[7:9], ["ϕ"; "ρ"; "σ"]) 
println("90% CI, Bayes")
prettyprint(bayes_inci[7:9], ["ϕ"; "ρ"; "σ"])

