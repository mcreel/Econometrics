function gmmresults()
    # example of GMM: draws from N(0,1). We estimate mean and variance.
    y = randn(1000,1)
    # 3 moment conditions implied by drawing from N(0,σ²):
    # mean = 0, variance = constant, skew = 0
    moments = θ -> [y.-θ[1] (y.^2.0).-θ[2] (y.-θ[1]).^3.0]
    # first round consistent
    W = Matrix{Float64}(I,3,3)
    θ = [0.0, 1.0]
    θhat, objvalue, D, ms, converged = gmm(moments, θ, W)
    # second round efficient
    Ωinv = inv(cov(ms))
    gmmresults(moments, θhat, Ωinv, "GMM example, two step")
    # CUE
    gmmresults(moments, θhat)
    return
end    

function gmmresults(moments, θ, weight="", title="", names="", efficient=true)
    n,g = size(moments(θ))
    CUE = false
    if weight !="" # if weight provided, use it
        θhat, objvalue, D, ms, converged = gmm(moments, θ, weight)
    else # do CUE
        θhat, objvalue, D, ms, weight, converged = gmm(moments, θ)
        CUE=true
    end
    k,g = size(D)
    # estimate variance
    efficient ? V = inv(D*weight*D')/n : V = V*D*weight*NeweyWest(ms)*weight*D'*V/n
    if names==""
        names = 1:k
        names = names'
    end
    se = sqrt.(diag(V))
    t = θhat ./ se
    p = 2.0 .- 2.0*cdf.(Ref(TDist(n-k)),abs.(t))
    PrintDivider()
    if title !="" printstyled(title, color=:yellow); println() end
    print("GMM Estimation Results    Convergence: ")
    printstyled(converged, color=:green)
    CUE ? message=" (CUE)" : message=" (user weight matrix)"
    printstyled(message, color=:cyan)
    println()
    println("Observations: ", n)
    println("Hansen-Sargan statistic: ", round(n*objvalue, digits=5))
    if g > k
        println("Hansen-Sargan p-value: ", round(1.0 - cdf(Chisq(g-k),n*objvalue), digits=5))
    end    
    a =[θhat se t p]
    println("")
    PrintEstimationResults(a, names)
    println()
    PrintDivider()
    return θhat, objvalue, V, D, weight, converged
end    
