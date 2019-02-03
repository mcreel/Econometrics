using Calculus
# moments should operate on a vector and return a nXg matrix
# of moment contributions

# ordinary GMM with provided weight matrix
# weight should be gXg positive definite
function gmm(moments, theta, weight)
    # average moments
    m = theta -> vec(mean(moments(theta),dims=1)) # 1Xg
    # moment contributions
    momentcontrib = theta -> moments(theta) # nXg
    # GMM criterion
    obj = theta -> ((m(theta))'weight*m(theta))
    # do minimization
    thetahat, objvalue, converged = fminunc(obj, theta)
    # derivative of average moments
    D = (Calculus.jacobian(m, vec(thetahat), :central))' 
    # moment contributions at estimate
    ms = momentcontrib(thetahat)
    return thetahat, objvalue, D, ms, converged
end

# CUE GMM, weight computed by NW
function gmm(moments, theta)
    # average moments
    m = theta -> vec(mean(moments(theta),dims=1)) # 1Xg
    # moment contributions
    momentcontrib = theta -> moments(theta) # nXg
    # weight
    weight = theta -> inv(NeweyWest(momentcontrib(theta)))
    # objective
    obj = theta -> m(theta)'*weight(theta)*m(theta)
    # do minimization
    thetahat, objvalue, converged = fminunc(obj, theta)
    # derivative of average moments
    D = (Calculus.jacobian(m, vec(thetahat), :central))' 
    # moment contributions at estimate
    ms = momentcontrib(thetahat)
    return thetahat, objvalue, D, ms, converged
end


# called with no args runs this example
function gmm()
    println("simple GMM example, sampling from N(0,1)")
    println("and using 3 moment conditions to estimate")
    println("do edit(gmm,()) to see the example code that this is running")
    # example of GMM: draws from N(0,1)
    y = randn(100,1)
    # 3 moment conditions
    moments = theta -> [y.-theta[1] (y.^2.0).-theta[2] (y.-theta[1]).^3.0]
    # first round consistent
    W = eye(3)
    theta = [0.0, 1.0]
    thetahat, objvalue, D, ms, converged = gmm(moments, theta, W)
    # second round efficient
    W = inv(cov(ms))
    thetahat, objvalue, D, ms, converged = gmm(moments, thetahat, W)
    # CUE
    thetahat, objvalue, D, ms, converged = gmm(moments, thetahat)
end    

