using Calculus, ForwardDiff

"""
    mle(model, θ, vc=1)

Maximum likelihood estimation.

Arguments:
- `model`: A function returning the vector of log-likelihoods for observations.
- `θ`: Initial parameter vector.
- `vc`: Method for computing variance; default is sandwich (1).

Returns:
- Estimated parameters, log-likelihood value, variance-covariance matrix, and convergence flag.
"""
function mle(model, θ, vc=1)
    avg_obj = θ -> -mean(vec(model(θ)))  # Negative average log-likelihood
    thetahat, objvalue, converged = fminunc(avg_obj, θ)  # Optimize

    objvalue = -objvalue  # Restore positive log-likelihood
    obj = θ -> vec(model(θ))  # Raw log-likelihood vector
    n = length(obj(θ))  # Number of observations
    k = length(θ)       # Number of parameters

    # Compute Jacobian (score contributions) using automatic diff or fallback
    scorecontrib = try
        ForwardDiff.jacobian(obj, thetahat)
    catch
        Calculus.jacobian(obj, thetahat, :central)
    end

    # Information matrix: I = Σ (∇ℓᵢ ⋅ ∇ℓᵢᵀ)
    I = scorecontrib' * scorecontrib

    # Hessian (expected or observed information matrix)
    H = try
        ForwardDiff.hessian(avg_obj, thetahat)
    catch
        Calculus.hessian(avg_obj, thetahat, :central)
    end

    # Compute variance-covariance matrix based on selected method
    vhat = if vc == 0
        inv(H)
    elseif vc == 1
        inv(H) * I * inv(H) / n
    elseif vc == 2
        inv(I)
    else
        error("Unknown vc option: $vc")
    end

    return thetahat, objvalue, vhat, converged
end



# with no args, runs example (also used for test)
function mle()
    x = [ones(10) 1:10]
    y = [1:2;3;4;5;6;7;8;9;10]
    β = zeros(2)
    model = β -> poisson(β, y, x)
    βhat, junk, junk, junk = mle(model, β)
    βhat
end
