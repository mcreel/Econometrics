using Statistics, Calculus, LinearAlgebra, DelimitedFiles

# compute moving average using p most recent values, including current value
function ma(x, p)
    m = similar(x)
    for i = p:size(x,1)
        m[i] = mean(x[i-p+1:i])
    end
    return m
end

# auxiliary model: HAR-RV
# Corsi, Fulvio. "A simple approximate long-memory model
# of realized volatility." Journal of Financial Econometrics 7,
# no. 2 (2009): 174-196.
function HAR(y)
    ylags = lags(y,10)
    X = [ones(size(y,1)) ylags[:,1]  mean(ylags[:,1:4],dims=2) mean(ylags[:,1:10],dims=2)]
    # drop missings
    y = y[11:end]
    X = X[11:end,:]
    βhat = X\y
    σhat = std(y-X*βhat)     
    vcat(βhat,σhat)
end

function aux_stat(y)
    y = abs.(y)
    m = mean(y)
    s = std(y)
    y = (y .- m)./s
    # look for evidence of volatility clusters
    mm = ma(y,5)
    mm = mm[5:end]
    clusters5 = quantile(mm,0.75)/quantile(mm, 0.25)
    mm = ma(y,10)
    mm = mm[10:end]
    clusters10 = quantile(mm,0.75)/quantile(mm, 0.25)
    ϕ = HAR(y)
    vcat(m, s, clusters5, clusters10, ϕ)
end

# map R^3 to param space: not used, but maybe will be
function θ2ϕ(θ)
    ϕ = similar(θ)
    ϕ[1] = exp(θ[1])
    ϕ[2] = 1.0 / (1.0 + exp(-θ[2]))
    ϕ[3] = exp(θ[3])
    return ϕ
end

# map param space to R^3: not used but maybe will be
function ϕ2θ(ϕ)
    θ = similar(ϕ)
    θ[1] = log(ϕ[1])
    θ[2] = -log(1.0/ϕ[2] - 1.0) 
    θ[3] = log(ϕ[3])
    return θ
end

# simple discrete time SV model
function SVmodel(n, θ, randdraws = 0.0, savedata=false)
    σe = θ[1] # mean of log variance: positive
    ρ = θ[2] # AR of log variance: between 0 and 1
    σu = θ[3] # std. dev. of log variance: positive
    burnin = 100
    ys = zeros(n)
    hlag = 0.0
    # generate shocks if not provided
    if randdraws == 0.0
        randdraws = randn(burnin+n,2)
    end
    @inbounds for t = 1:burnin+n
        h = ρ*hlag + σu*randdraws[t,1] # log variance follows AR(1)
        y = σe*exp(h./2.0)*randdraws[t,2]
        if t > burnin
            ys[t-burnin] = y
        end
        hlag = h
    end
    if savedata
        writedlm("svdata.txt", ys)
    end    
    sqrt(n)*aux_stat(ys)
end

function SVmoments(m, n, θ, randdraws)
    S = size(randdraws, 1)
    mm = zeros(S,size(m,1))
    for s=1:S
        mm[s,:] = SVmodel(n, θ, randdraws[s,:,:]) - m
    end
    return mm
end    

function main()
# generate the sample
θtrue = [exp(-0.736/2.0), 0.9, 0.363] # true param values, on param space
lb = [0.0, 0.0, 0.0]
ub = [3.0, 1.0, 2.0]
n = 500 # sample size
burnin = 100
randdraws = randn(n+burnin,2)
m0 =  SVmodel(n, θtrue, randdraws, true) # generate the sample and save the data
# Estimation by indirect inference
S = 100 # number of simulation reps
randdraws = randn(S,n+burnin,2) # fix the shocks to control "chatter" (includes the burnin period)
ms = θ -> SVmoments(m0, n, θ, randdraws)
m = θ -> vec(mean(ms(θ),dims=1)) # 1Xg
weight = θ -> inv(cov(ms(θ)))
obj = θ -> m(θ)'weight(θ)*m(θ)
thetahat, objvalue, converged, details = samin(obj, θtrue, lb, ub; ns = 5, verbosity = 2, rt = 0.5)
# compute the estimated standard errors and CIs
D = (Calculus.jacobian(m, vec(thetahat), :central))
W = weight(thetahat)
V = inv(D'*W*D)
se = sqrt.(diag(V))

println("true values, estimates, st. error, and limits of 95% CI")
prettyprint([θtrue thetahat se thetahat-1.96*se thetahat+1.96*se],["true value", "estimate", "std. err.", "CI lower", "CI upper"])
end
main()
