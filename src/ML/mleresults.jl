"""
    mleresults(model, θ)

Do maximum likelihood estimation, where model is a function that gives
the vector of log likelihoods of the observations, and θ is the
parameter vector. call edit(mleresults()) to see a coded example.

"""
function mleresults()
    println("execute edit(mleresults,()) to examine the example code")
    x = [ones(10) 1:10]
    y = [1:2;3;4;5;6;7;8;9;10]
    β = zeros(2)
    model = β -> poisson(β, y, x)
    βhat, junk, junk, junk = mleresults(model, β, "simple ML example")
    βhat
end

function mleresults(model, θ, title="", names=""; vc=1)
    n = size(model(θ),1)
    thetahat, objvalue, V, converged = mle(model, θ, vc)
    k = size(V,1)
    if names==""
        names = 1:k
        names = names'
    end
    se = sqrt.(diag(V))
    t = thetahat ./ se
    p = 2.0 .- 2.0*cdf.(Ref(TDist(n-k)),abs.(t))
    PrintDivider()
    if title !="" printstyled(title, color=:yellow); println() end
    print("MLE Estimation Results    Convergence: ")
    printstyled(converged, color=:green)
    println()
    println("Average Log-L: ", round(objvalue; digits=5), "   Observations: ", n)
    if vc==1
        println("Sandwich form covariance estimator")
    elseif vc==2
        println("Inverse Hessian form covariance estimator")
    else
        println("OPG form covariance estimator")
    end
    a =[thetahat se t p]
    println("")
    PrintEstimationResults(a, names)
    println()
    println("Information Criteria")
    caic = -2.0*n*objvalue + k*(log(n)+1.0)
    bic = -2.0*n*objvalue + k*log(n)
    aic = -2.0*n*objvalue + 2.0*k
    infocrit = [caic; bic; aic]
    infocrit = round.([infocrit infocrit/n],digits=5)
    clabels = ["Crit.", "Crit/n"]
    rlabels = ["CAIC ", "BIC ", "AIC "]
    prettyprint(infocrit, clabels, rlabels)
    PrintDivider()
    return thetahat, objvalue, V, converged
end
