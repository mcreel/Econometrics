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
    x = [ones(10) 1:10]
    y = [1:2;3;4;5;6;7;8;9;10]
    β = zeros(2)
    # 3 moment conditions
    moments = β -> x.*(y - exp.(x*β))
    # first round consistent
    W = eye(2)
    βhat, objvalue, D, ms, converged = gmm(moments, β, W)
    βhat
end    

