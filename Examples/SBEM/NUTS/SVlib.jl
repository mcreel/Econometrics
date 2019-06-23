using Statistics, Random, Distributions

# returns the variable (or matrix), lagged p times,
# with the first p rows filled with ones (to avoid divide errors)
# remember to drop those rows before doing analysis
function lag(x,p::Int64)
	n = size(x,1)
	lagged_x = zeros(eltype(x),n,p)
	lagged_x = [ones(p); x[1:n-p]]
end

# lags of a vector from 1 to p in a p column array
function  lags(x,p)
	n = size(x,1)
	lagged_x = zeros(eltype(x),n,p)
	for i = 1:p
		lagged_x[:,i] = lag(x,i)
	end
    return lagged_x
end

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
    y, m, s  = stnorm(abs.(y))
    # look for evidence of volatility clusters
    mm = ma(y,5)
    mm = mm[5:end]
    clusters = quantile(mm,0.75) -quantile(mm, 0.25)
    ϕ = HAR(y)
    vcat(m, s, clusters, ϕ)
end

# the dgp: simple discrete time stochastic volatility (SV) model
function SVmodel(σe, ρ, σu, n, shocks_u, shocks_e)
    burnin = size(shocks_u,1) - n
    hlag = 0.0
    h = ρ.*hlag .+ σu.*shocks_u[1] # figure out type
    y = σe.*exp(h./2.0).*shocks_e[1]
    ys = zeros(eltype(y),n)
    for t = 1:burnin+n
        h = ρ.*hlag .+ σu.*shocks_u[t]
        y = σe.*exp(h./2.0).*shocks_e[t]
        if t > burnin 
            ys[t-burnin] = y
        end    
        hlag = h
    end
    #plot(ys)
    sqrt(n)*aux_stat(ys)
end

