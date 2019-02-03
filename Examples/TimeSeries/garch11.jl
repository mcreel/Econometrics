# computes GARCH(1,1) log likelihood
using Statistics
function garch11(theta, y)
    # dissect the parameter vector
    mu = theta[1]
    rho = theta[2]
    omega = theta[3]
    alpha = theta[4]
    beta = theta[5]
    resid = y[2:end] .- mu .- rho*y[1:end-1]
    n = size(resid,1)
    h = zeros(n)
    # initialize variance; either of these next two are reasonable choices
    h[1] = var(y[1:10])
    #h[1] = var(y)
    rsq = resid.^2.0
    for t = 2:n
        h[t] = omega + alpha*rsq[t-1] + beta*h[t-1]
    end
    logL = -log(sqrt(2.0*pi)) .- 0.5*log.(h) .- 0.5*rsq./h
end    
