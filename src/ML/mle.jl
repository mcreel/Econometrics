using Calculus, ForwardDiff
"""
    mle(model, θ, vc=" ")

Maximum likelihood estimation.

* model should be a function that gives the vector of log likelihoods of
the observations
* θ is the parameter vector
* vc is the method for computing the variance, default is sandwich
* execute mleresults() for an example, or edit(mleresults,()) to see the code.

"""

function mle(model, θ, vc=1)
    avg_obj = θ -> -mean(vec(model(θ))) # average log likelihood
    thetahat, objvalue, converged = fminunc(avg_obj, θ) # do the minimization of -logL
    objvalue = -objvalue
    obj = θ -> vec(model(θ)) # unaveraged log likelihood
    n = size(obj(θ),1) # how many observations?
    # derivatives, use automatic, with numeric as fallback
    # information matrix
    scorecontrib = try
        ForwardDiff.jacobian(obj, vec(thetahat))
    catch
        Calculus.jacobian(obj, vec(thetahat), :central)
    end   
    k = size(θ,1)
    I = zeros(k,k)
    for i = 1:n
    	I .+= scorecontrib[i,:]*scorecontrib[i,:]'
    end
    I ./= n
    # Hessian
    J = try
        ForwardDiff.hessian(avg_obj, vec(thetahat))
    catch
        Calculus.hessian(avg_obj, vec(thetahat), :central)
    end
    Jinv = inv(J)
    # three ways of computing covariance estimates
    if vc==2
        V = Jinv/n      # other possibilities
    elseif vc==3
        V = inv(I)/n
    else
        V = Jinv*I*Jinv/n # sandwich form is preferred, and is the default
    end
    return thetahat, objvalue, V, converged
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
