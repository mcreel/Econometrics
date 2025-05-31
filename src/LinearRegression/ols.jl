using StatsFuns, StatsModels, DataFrames, LinearAlgebra

"""
    ols(y,x)

Compute the ordinary least squares regresssion coefficients, and report results. Can also compute OLS subject to linear restrictions.

ols() with no arguments will run an example, execute edit(ols,()) to see the code.

"""
function ols()
    y = rand(100,1)
    x = [ones(100,1) rand(100,3)]
    names = ["const" "x2"  "x3" "x4"]
    ols(y,x)
    ols(y,x,names=names)
    ols(y,x,names=names, vc="nw")
    b = ols(y,x,silent=true)[1]
    show(b)
    # restricted LS
    R = [1 0 0 0]
    r = 1
    ols(y,x,R=R,r=r)
    # using dataframe
    println("using a formula and dataframe")
    d = DataFrame(x=rand(10),y=rand(10),z=rand(10))
    f = @formula(y ~ 1+x+z+x*z)
    ols(f,d)
    return
end


# version for dataframe and formula
function ols(f::FormulaTerm, df::DataFrame; R=[], r=[], names="", vc="white", silent=false)
    y = modelmatrix(f.lhs,df)
    x = modelmatrix(f,df)
    ff = apply_schema(f,schema(df))
    ols(y, x, names=coefnames(ff.rhs))
end

function ols(y::Array{Float64}, x::Array{Float64,2}; R=[], r=[], names="", vc="white", silent=false)
    n,k = size(x)
    if names==""
        names = 1:k
        names = names'
    end
    b, fit, e = lsfit(y,x)
    df = n-k
    sigsq = (e'*e/df)[1,1]
    xx_inv = inv(x'*x)
    ess = (e' * e)[1,1]
    # Restricted LS?
    if R !=[]
        q = size(R,1)
        P_inv = inv(R*xx_inv*R')
        b = b .- xx_inv*R'*P_inv*(R*b.-r)
        e = y-x*b;
        ess = (e' * e)[1,1]
        df = n-k-q
        sigsq = ess/df
        A = Matrix{Float64}(I, k, k) .- xx_inv*R'*P_inv*R;  # the matrix relating b and b_r
    end
    # Ordinary or het. consistent variance estimate
    if vc=="white"
        xe = x.*e
        varb = xx_inv*xe'xe*xx_inv
    elseif vc=="nw"
        xe = x.*e
        lags = Int(max(round(n^0.25),1.0))
        varb = n*xx_inv*NeweyWest(xe,lags)*xx_inv
    else
        varb= xx_inv.*sigsq
    end
    # restricted LS?
    if R !=[]
        varb = A*varb*A'
    end
    # common to both ordinary and restricted
    seb = sqrt.(diag(varb))
    seb = seb.*(seb.>1e-16) # round off to zero when there are restrictions
    t = b ./ seb
    tss = y .- mean(y)
    tss = (tss'*tss)[1,1]
    rsq = (1.0 - ess / tss)
    if !silent
        println()
        PrintDivider()
        if R==[]
            @printf("  OLS estimation, %d observations\n", n)
        else
            @printf("  Restricted LS estimation, %d observations\n", n)
        end
        @printf("  R²: %f   σ²: %f\n", rsq, sigsq)
        p = 2.0 .- 2.0*tdistcdf.(df, abs.(t))
        results = [b seb t p]
        if vc=="white"
            println("  White's covariance estimator")

        elseif vc=="nw"
            println("  Newey-West covariance estimator")
        end
        println()
        PrintEstimationResults(results, names)
        PrintDivider()
    end
    return b, varb, e, ess, rsq
end
