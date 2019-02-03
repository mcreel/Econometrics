using StatsFuns, LinearAlgebra

"""
    tsls(y, x, z)

    Two stage least squares regression of y on x, using instruments z

"""
function tsls(y, x, z; names="", vc="white", silent=false)
    n,k = size(x)
    if names==""
        names = 1:k
        names = names'
    end
    xhat = z*(z\x)
    xx_inv = inv(xhat'x)
    b = xx_inv*xhat'y
    e = y - x*b
    ess = (e'e)[1,1]
    sigsq = ess/(n-k)

    # Ordinary or het. consistent variance estimate
    if vc=="white"
        xe = xhat.*e
        varb = xx_inv*xe'xe*xx_inv
    elseif vc=="nw"
        xe = xhat.*e
        lags = Int(max(round(n^0.25),1.0))
        varb = n*xx_inv*NeweyWest(xe,lags)*xx_inv
    else
        varb= xx_inv.*sigsq
    end
    seb = sqrt.(diag(varb))
    t = b ./ seb
    tss = y .- mean(y)
    tss = (tss'*tss)[1,1]
    rsq = (1.0 - ess / tss)
    labels = ["coef", "se", "t", "p"]
    if !silent
        PrintDivider()
        @printf("  2SLS estimation, %d observations\n", n)
        @printf("  R²: %f   σ²: %f\n", rsq, sigsq)
        p = 2.0 .- 2.0*tdistcdf.(n-k, abs.(t))
        results = [b seb t p]
        if vc=="white"
            println("  White's covariance estimator")

        elseif vc=="nw"
            println("  Newey-West covariance estimator")
        end
        println()
        prettyprint(results, labels, names)
        PrintDivider()
    end
    return b, varb, e, ess, rsq
end
