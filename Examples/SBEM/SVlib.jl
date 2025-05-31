# The library of functions for the simple SV model

using Statistics, Random

# method that generates the sample
function auxstat(θ, reps)
    stats = zeros(reps,11)
    if InSupport(θ)
        n = 500
        burnin = 100
        for rep = 1:reps
            y = SVmodel(θ, n, burnin)
            stats[rep,:] = auxstat(y)
        end
    end    
    stats
end

# method that uses pre-generated shocks
function auxstat(θ, ϵ, u)
    ϕ = θ[1]
    ρ = θ[2]
    σ = θ[3]
    reps = size(ϵ,2)
    stats = zeros(reps,11)
    if InSupport(θ)
        n = 500
        burnin = 100
        for rep = 1:reps
            hlag = 0.0
            ys = zeros(n)
            @inbounds for t = 1:burnin+n
                h = ρ*hlag + σ*u[t, rep]
                y = ϕ*exp(h/2.0)*ϵ[t, rep]
                if t > burnin 
                    ys[t-burnin] = y
                end    
                hlag = h
            end
            stats[rep,:] = auxstat(ys)
        end
    end    
    stats
end

# method for a given sample
function auxstat(y)
	s = std(y)
	y = abs.(y)
	m = mean(y)
	s2 = std(y)
	y = y ./ s2
	k = std((y).^2.0)
	c = cor(y[1:end-1],y[2:end])
	# ratios of quantiles of moving averages to detect clustering
	q = try
	    q = quantile((ma(y,3)[3:end]), [0.25, 0.75])
	catch
	    q = [1.0, 1.0]
	end
	c1 = log(q[2]/q[1])
	stats = sqrt(size(y,1)) .* vcat(m, s, s2, k, c, c1, HAR(y))
end

function SVmodel(θ, n, burnin)
    ϕ = θ[1]
    ρ = θ[2]
    σ = θ[3]
    hlag = 0.0
    ys = zeros(n)
    if InSupport(θ)
        @inbounds for t = 1:burnin+n
            h = ρ*hlag + σ*randn()
            y = ϕ*exp(h/2.0)*randn()
            if t > burnin 
                ys[t-burnin] = y
            end    
            hlag = h
        end
    end    
    ys
end

function TrueParameters()
    [exp(-0.736/2.0), 0.9, 0.363]
end

function PriorSupport()
    lb = [0.05, 0.0, 0.05]
    ub = [2.0, 0.999, 1.0]
    lb,ub
end    

# prior checks that we're in the bounds, and that the unconditional std. dev. of log vol is not too high
# returns 1 if this is true, zero otherwise. Value is not important, as it's constant
function Prior(θ)
    InSupport(θ) ? 1.0 : 0.0
end

# check if parameter is in support. In this case, we require
# the bounds, and that the unconditional variance of the volatility
# shock be limited
function InSupport(θ)
    lb,ub = PriorSupport()
    if all(θ .>= lb) & all(θ .<= ub)
	return  (θ[3]/sqrt(1.0 - θ[2]^2.0) < 5.0)
    else
    return false
    end	
end

function PriorDraw()
    lb, ub = PriorSupport()
    ok = false
    θ = 0.0
    while !ok
        θ = (ub-lb).*rand(size(lb,1)) + lb
        ok = InSupport(θ)
    end
    return θ
end    

# taken from https://github.com/mcreel/Econometrics  
# # returns the variable (or matrix), lagged p times,
# with the first p rows filled with ones (to avoid divide errors)
# remember to drop those rows before doing analysis
function lag(x::Array{Float64,2},p::Int64)
	n,k = size(x)
	lagged_x = [ones(p,k); x[1:n-p,:]]
end

function lag(x::Array{Float64,1},p::Int64)
	n = size(x,1)
	lagged_x = [ones(p); x[1:n-p]]
end	

# taken from https://github.com/mcreel/Econometrics  
# returns the variable (or matrix), lagged from 1 to p times,
# with the first p rows filled with ones (to avoid divide errors)
# remember to drop those rows before doing analysis
function  lags(x::Array{Float64,2},p)
	n, k = size(x)
	lagged_x = zeros(eltype(x),n,p*k)
	for i = 1:p
		lagged_x[:,i*k-k+1:i*k] = lag(x,i)
	end
    return lagged_x
end	

function  lags(x::Array{Float64,1},p)
	n = size(x,1)
	lagged_x = zeros(eltype(x), n,p)
	for i = 1:p
		lagged_x[:,i] = lag(x,i)
	end
    return lagged_x
end

# taken from https://github.com/mcreel/Econometrics  
# compute moving average using p most recent values, including current value
function ma(x, p)
    m = zeros(size(x))
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
