#= note: to run this in parallel, from the system
prompt, and before starting julia, do
export JULIA_NUM_THREADS=X
where X is the desired number of threads. On
my 32 core system, 10 is good (cache contention issues)
Remember to re-set threads to 1 before using MPI
=#
using Statistics, Econometrics, LinearAlgebra, Plots
include("QIVmodel.jl")
function main()
LNW, X, Z = getdata()
n = size(LNW,1)
mcmcreps = 100000
burnin = 10000
# use st. errs. from ordinary QR as quide to setting tuning.
basetuning = [0.1, 0.05, 0.03, 0.005, 0.3, 0.2, 0.2]
# adjust this up or down to achieve decent acceptance rate
scale = 0.1 # higher for higher taus 
tuning = scale*basetuning
EducEffect = zeros(9,3)
Constant = zeros(9,3)
for i = 1:9
    τ = round(i/10,digits=1)
    Σ = τ*(1.0-τ)*Z'Z/n
    Σinv = inv(Σ)
    θ = X\LNW  # OLS start values
    # to do ordinary QR via MCMC, set Z=X (just to verify that MCMC works!)
    lnL = θ -> likelihood(θ, LNW, X, Z, τ, Σinv)
    Proposal = θ -> proposal(θ, tuning)
    Prior = θ -> prior(θ)
    chain = mcmc(θ, mcmcreps, burnin, Prior, lnL, Proposal)
    V = cov(chain[:,1:end-1])
    d = dstats(chain)
    Constant[i,:] = d[1,[1,7,8]] # save posterior mean, 5% and 95% quantile
    EducEffect[i,:] = d[2,[1,7,8]] # save posterior mean, 5% and 95% quantile
    #scatter(chain[:,1],chain[:,2], show=true, reuse=false)
end

# plot results
# educ
τ = (1:9)/10
plot(reuse=false)
educ = EducEffect[:,1]
lb = EducEffect[:,2]
ub = EducEffect[:,3]
plot!(τ, [educ educ], fillrange=[lb ub], fillalpha=0.3, c=:green, legend=:none)
xticks!((1:9)/10)
#savefig("EducJUNK.svg")
# constant
τ = (1:9)/10
plot(reuse=false)
c = Constant[:,1]
lb = Constant[:,2]
ub = Constant[:,3]
plot!(τ, [c c], fillrange=[lb ub], fillalpha=0.3, c=:green, legend=:none)
xticks!((1:9)/10)
#savefig("ConstantJUNK.svg")
return
end
main()
