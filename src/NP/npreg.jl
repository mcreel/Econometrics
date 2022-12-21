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

using Random, Econometrics, Plots
function npreg(bandwidth=-1, silent=false)
    itworked = true
    try
        if !silent
            println("npreg(), with no arguments, runs a simple example")
            println("npreg(bw) will run the example with your chosen bandwidth (Float64)")
            println("execute edit(npreg,()) to see the code")
        end
        k = 1 # number of regressors
        Random.seed!(1) # set seed to 
        n = 1000
        bandwidth == -1 ? bandwidth = 0.25*n^(-1.0/(4 + k)) : nothing
        neval = 100
        x = rand(n)*pi*2.0 # from 0 to 2π   
        xeval = collect(range(pi/2., stop=pi*1.5, length=neval)) # from π/2 to 3π/2  
        y = cos.(x) .+ 0.25 .* cos.(3.0.*x) + 0.1*randn(n)
        # npreg wants args to be matrices, not vectors, so convert them
        y = reshape(y,n,1)
        x = reshape(x,n,1)
        xeval = reshape(xeval,neval,1)
        ytrue = cos.(xeval) .+ 0.25.*cos.(3.0*xeval)
        weights = kernelweights(x, xeval, bandwidth, true, "knngaussian", 200)
        yhat, y50, y05, y95 = npreg(y, x, xeval, weights, order=1, do_median=true, do_ci=true)
        labels = ["true" "mean" "median" "0.05 quantile" "0.95 quantile"]
        title = "Kernel regression and quantiles"
        p = Plots.plot(xeval, [ytrue yhat y50 y05 y95], labels=labels, title=title)
        display(p)
    catch
        itworked = false
    end
    itworked # report this for testing
end 

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
