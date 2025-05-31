using ForwardDiff, Calculus
# moments should operate on a vector and return a nXg matrix
# of moment contributions

# ordinary GMM with provided weight matrix
# weight should be gXg positive definite
function gmm(moments, θ, weight)
    # average moments
    m = θ -> vec(mean(moments(θ),dims=1)) # 1Xg
    # moment contributions
    # GMM criterion
    obj = θ -> dot(m(θ), weight, m(θ))
    # do minimization
    θhat, objvalue, converged = fminunc(obj, θ)
    # derivative of average moments
    D = try
        ForwardDiff.jacobian(m, vec(θhat))'  # jacobian get ∂m/∂θ', but D is the transpose
    catch    
        Calculus.jacobian(m, vec(θhat), :central)'
    end    
    # moment contributions at estimate
    ms = moments(θhat)
    return θhat, objvalue, D, ms, converged
end

# CUE GMM, weight computed by NW
function gmm(moments, θ)
    # average moments
    m = θ -> vec(mean(moments(θ),dims=1)) # 1Xg
    # weight
    Ωinv = θ -> inv(NeweyWest(moments(θ)))
    # objective
    obj = θ -> dot(m(θ),Ωinv(θ),m(θ))
    # do minimization
    θhat, objvalue, converged = fminunc(obj, θ)
    # derivative of average moments
    D = try
        ForwardDiff.jacobian(m, vec(θhat))' 
    catch
        Calculus.jacobian(m, vec(θhat), :central)' 
    end
    # moment contributions at estimate
    ms = moments(θhat)
    return θhat, objvalue, D, ms, Ωinv(θhat), converged
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

