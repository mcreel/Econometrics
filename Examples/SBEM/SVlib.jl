using Econometrics, Statistics, Random, Distributions, DelimitedFiles

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
    y = abs.(y) # add 1 to avoid log(0)
    m = mean(y)
    s = std(y)
    # look for evidence of volatility clusters
    mm = ma(y,5)
    mm = mm[5:end]
    clusters5 = 0.0
    clusters10 = 0.0
    try
        clusters5 = quantile(mm,0.75)/quantile(mm, 0.25)
    catch
        clusters5 = 1.0
    end    
    mm = ma(y,10)
    mm = mm[10:end]
    try
        clusters10 = quantile(mm,0.75)/quantile(mm, 0.25)
    catch
        clusters10 = 1.0
    end    
    ϕ = HAR(y)
    vcat(m, s, clusters5, clusters10, ϕ)
end

# the dgp: simple discrete time stochastic volatility (SV) model
function SVmodel(σe, ρ, σu, n, shocks_u, shocks_e, savedata=false)
    burnin = size(shocks_u,1) - n
    hlag = 0.0
    h = ρ.*hlag .+ σu.*shocks_u[1] # figure out type
    y = σe.*exp(h./2.0).*shocks_e[1]
    ys = zeros(n,1)
    for t = 1:burnin+n
        h = ρ.*hlag .+ σu.*shocks_u[t]
        y = σe.*exp(h./2.0).*shocks_e[t]
        if t > burnin 
            ys[t-burnin] = y
        end    
        hlag = h
    end
    if savedata == true
        writedlm("svdata.txt", ys)
    end    
    sqrt(n)*aux_stat(ys)
end

# for GMM 
function SVmoments(m, n, θ, randdraws)
    S = size(randdraws, 1)
    mm = zeros(S,size(m,1))
    σe = θ[1]
    ρ = θ[2]
    σu = θ[3]
    for s=1:S
        mm[s,:] = SVmodel(σe, ρ, σu, n, randdraws[s,:,1], randdraws[s,:,2])
    end
    mm .- m' # both are already scaled by root-n
end

# asymptotic Gaussian likelihood function of statistic
function logL(θ, m, n, shocks_u, shocks_e)
    σe = θ[1]
    ρ = θ[2]
    σu = θ[3]
    S = size(shocks_u,2)
    k = size(m,1)
    ms = zeros(eltype(SVmodel(σe, ρ, σu, n, shocks_u[:,1], shocks_e[:,1])), S, k)
    # this loop could be parallelized!
    for s = 1:S
        ms[s,:] = SVmodel(σe, ρ, σu, n, shocks_u[:,s], shocks_e[:,s])
    end
    mbar = mean(ms,dims=1)[:]
    Σ = cov(ms)
    x = (m .- mbar)
    logL = try
        logL = -0.5*log(det(Σ)) - 0.5*x'*inv(Σ)*x
    catch
        logL = -Inf
    end
end

# uniform random walk, with bounds check
function proposal1(current, tuning, lb, ub)
    trial = copy(current)
    if rand() > 0.1
        i = rand(1:size(current,1))
        trial[i] = current[i] + tuning[i].*randn()
    else
        trial = lb + (ub - lb).*rand(size(lb))
    end    
    return trial
end

function proposal2(current, cholV, lb, ub)
    trial = copy(current)
    if rand() > 0.1
        trial += cholV'*randn(size(trial))
    else
        trial = lb + (ub - lb).*rand(size(lb))
    end


end

function prior(theta, lb, ub)
    a = 0.0
    if(all((theta .>= lb) .& (theta .<= ub)))
        a = 1.0
    end
    return a
end

