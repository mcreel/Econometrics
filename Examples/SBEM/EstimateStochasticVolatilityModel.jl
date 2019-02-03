using Statistics, StatsPlots
pyplot(show=true, reuse=false)
# map R^3 to param space
function θ2ϕ(θ)
    ϕ = similar(θ)
    ϕ[1] = exp(θ[1])
    ϕ[2] = 1.0 / (1.0 + exp(-θ[2]))
    ϕ[3] = exp(θ[3])
    return ϕ
end
# map param space to R^3
function ϕ2θ(ϕ)
    θ = similar(ϕ)
    θ[1] = log(ϕ[1])
    θ[2] = -log(1.0/ϕ[2] - 1.0) 
    θ[3] = log(ϕ[3])
    return θ
end
# simple discrete time SV model
function SVmodel(n, θ, randdraws = 0.0)
    ϕ = θ2ϕ(θ)
    α = ϕ[1] # mean of log variance: positive
    ρ = ϕ[2] # AR of log variance: between 0 and 1
    σᵤ = ϕ[3] # std. dev. of log variance: positive
    burnin = 1000
    ys = zeros(n)
    hlag = 0.0
    # generate shocks if not provided
    if randdraws == 0.0
        randdraws = randn(burnin+n,2)
    end
    @inbounds for t = 1:burnin+n
        h = α + ρ*(hlag-α) + σᵤ*randdraws[t,1] # log variance follows AR(1)
        y = exp(h/2.0)*randdraws[t,2]
        if t > burnin
            ys[t-burnin] = y
        end
        hlag = h
    end
    σ = exp(hlag/2.0) # the latent volatility in last period
    return ys, σ
end

# auxiliary model: HAR-RV(p)
# Corsi, Fulvio. "A simple approximate long-memory model
# of realized volatility." Journal of Financial Econometrics 7,
# no. 2 (2009): 174-196.
function HAR(y,p)
    σhat = std(y)
    y ./= σhat
    RV = abs.(y)
    RVlags = lags(RV,p)
    X = [ones(size(y,1)) RVlags]
    # drop missings
    RV = RV[p+1:end]
    X = X[p+1:end,:]
    βhat = X\RV
    return vcat(βhat,σhat)
end

function II_moments(θ, ϕhat, n, randdraws)
    S = size(randdraws,1)
    n_auxparams = size(ϕhat,1)
    p = n_auxparams - 2 # number of lags in HAR
    ϕhatS = zeros(S, n_auxparams)
    for s = 1:S
        yₛ, junk = SVmodel(n, θ, randdraws[s,:,:])
        ϕhatS[s,:] = HAR(yₛ,p)
    end
    ms = ϕhat' .- ϕhatS
    return  ms # the moments, in a SxG matrix
end

function main()
# generate the sample
ϕ₀ = [0.7, 0.95, 0.2] # true param values, on param space
θₒ  = ϕ2θ(ϕ₀) # true parameter values, on R^2
n = 1000 # sample size
y, σ = SVmodel(n, θₒ) # generate the sample
plot(y)
density(y)

# Estimation by indirect inference
S = 20 # number of simulation reps
randdraws = randn(S,n+1000,2) # fix the shocks to control "chatter" (includes the burnin period)
p = 4
ϕhat = HAR(y,p) # statistics using real sample
moments = θ -> II_moments(θ, ϕhat, n, randdraws)
thetahat, junk, junk, junk, junk = gmm(moments, θₒ, 1.0)
gmmresults(moments, thetahat, "", "This minimizes m(θ)'inv[Ω(θ)]m(θ), but std errs, etc. are not correct")
θ2ϕ(thetahat)
end
main()
