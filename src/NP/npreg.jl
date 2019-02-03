using Random
"""
nonparametric regression, using pre-generated weights
can be local constant, local linear, or local quadratic

y is n X p

x is n X dimx

xeval is neval X dimx

weights is n X neval (different weights for each eval. point)

weights should sum to one by columns, they are typical 
nonparametric weights from a kernel.

execute npreg() for an example
"""
function npreg()
    println("npreg(), with no arguments, runs a simple example")
    println("execute edit(npreg,()) to see the code")
    #using Plots
    order = 2 # 0 for local constant, 1 for local linear, 2 for local quadratic
    k = 3 # number of regressors
    Random.seed!(1) # set seed to enable testing
    bandwidth = 0.07
    n = 100000
    neval = 100
    x = rand(n,k)*pi*2.0
    xeval = [pi*ones(neval,k-1) range(pi/2., stop=pi*1.5, length=neval)]
    y = cos.(sum(x,dims=2)) + cos.(2. * sum(x,dims=2)) + 0.5*randn(n,1)
    ytrue = cos.(sum(xeval,dims=2)) + cos.(2. * sum(xeval,dims=2))
    weights = kernelweights(x, xeval, bandwidth, true, "gaussian", 200)
    yhat = npreg(y, x, xeval, weights, order)
    #plot(xeval[:,k], [yhat ytrue])
end 

function npreg(y::Array{Float64,2}, x::Array{Float64,2}, xeval::Array{Float64,2}, weights::Array{Float64,2}, order::Int64=1, do_median::Bool=false, do_ci::Bool=false)
    weights = sqrt.(weights)
    neval, dimx = size(xeval)
    n, dimy = size(y)
    # local constant
    if order==0            
        X = ones(n,1)
        Xeval = ones(neval,1)
    elseif order==1    
    # local linear    
        X = [ones(n,1) x]
        Xeval = [ones(neval,1) xeval]
    else
    # local quadratic
        stacked = [x; xeval]
        # cross products
        CP = zeros(n+neval, Int((dimx-1)*dimx/2))
        cpind = 0
        for i = 1:dimx-1
            for j = (i+1):dimx
                cpind += 1
                CP[:,cpind] = stacked[:,[i]].*stacked[:,[j]]
            end
        end
        ZZ = [ones(n+neval,1) stacked CP]
        X = view(ZZ,1:n,:)
        Xeval = view(ZZ,(n+1):n+neval,:)
    end
    # do the fit
    yhat = zeros(neval, dimy)
    for i = 1:neval
        WX = weights[:,i] .* X
        Wy = weights[:,i] .* y
        b = WX\Wy
        yhat[i,:] = Xeval[[i],:]*b
    end    
    return yhat
end
