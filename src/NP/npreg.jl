using Random, Econometrics
"""
nonparametric regression and quantile regression
* uses pre-generated kernel weights
* can be local constant, local linear, or local quadratic

y is n X p

x is n X dimx

xeval is neval X dimx

weights is n X neval (different weights for each eval. point)

weights should sum to one by columns, they are typical 
nonparametric weights from a kernel.

keyword order=0,1 or 2 for local consant, local linear, or local quadratic
(default local linear)

keyword do_median (default false)
keyword do_ci to compute .05 and 0.95 quantiles

execute npreg() for an example
"""
function npreg()
    println("npreg(), with no arguments, runs a simple example")
    println("execute edit(npreg,()) to see the code")
    #using Plots
    k = 3 # number of regressors
    Random.seed!(1) # set seed to enable testing
    n = 10000
    bandwidth = 0.25*n^(-1.0/(4 + k))
    neval = 100
    x = rand(n,k)*pi*2.0
    xeval = [pi*ones(neval,k-1) range(pi/2., stop=pi*1.5, length=neval)]
    y = cos.(sum(x,dims=2))
    y = y + 0.1*randn(size(y))
    ytrue = cos.(sum(xeval,dims=2))
    weights = kernelweights(x, xeval, bandwidth, true, "knngaussian", 200)
    yhat, y50, y05, y95 = npreg(y, x, xeval, weights, order=1, do_median=true, do_ci=true)
    println("true, mean, median, q05, q95")
    prettyprint([ytrue yhat y50 y05 y95])
    return yhat[1,1] # for testing
end 

#=
# to do testing, should use nonrandom data
function npreg(x::Int64)
    #println("npreg(), with a single integer argument, runs a test, using the npreg() example")
    #using Plots
    k = 3 # number of regressors
    Random.seed!(1) # set seed to enable testing
    n = 10000
    bandwidth = 0.25*n^(-1.0/(4 + k))
    neval = 100
    x = rand(n,k)*pi*2.0
    xeval = [pi*ones(neval,k-1) range(pi/2., stop=pi*1.5, length=neval)]
    y = cos.(sum(x,dims=2))
    y = y + 0.1*randn(size(y))
    ytrue = cos.(sum(xeval,dims=2))
    weights = kernelweights(x, xeval, bandwidth, true, "knngaussian", 200)
    yhat, y50, y05, y95 = npreg(y, x, xeval, weights, order=1, do_median=true, do_ci=true)
    return yhat[1,1] # for testing
end 
=#


function npreg(y, x, xeval, weights; order=1, do_median=false, do_ci=false)
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
    y50 = nothing
    y05 = nothing
    y95 = nothing
    if do_median y50 = zeros(neval, dimy) end
    if do_ci
        y05 = zeros(neval, dimy)
        y95 = zeros(neval, dimy)
    end
    for i = 1:neval
        WX = weights[:,i] .* X
        Wy = weights[:,i] .* y
        b = WX\Wy
        yhat[i,:] = Xeval[i:i,:]*b
        if do_median
            for j = 1:dimy
                y50[i,j] = (Xeval[i:i,:]*qreg_coef(Wy[:,j], WX, 0.5))[1,1]
            end
        end    
        if do_ci
            for j = 1:dimy
                y05[i,j] = (Xeval[i:i,:]*qreg_coef(Wy[:,j], WX, 0.05))[1,1]
                y95[i,j] = (Xeval[i:i,:]*qreg_coef(Wy[:,j], WX, 0.95))[1,1]
            end
        end    
    end
    if do_median==false
        return yhat
    elseif do_ci==false
        return yhat, y50
    else
        return yhat, y50, y05, y95
    end    
end
