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
    scorecontrib = Calculus.jacobian(obj, vec(thetahat), :central)
    I = cov(scorecontrib)
    J = Calculus.hessian(avg_obj, vec(thetahat), :central)
    Jinv = inv(J)
    if vc==2
        V = Jinv/n      # other possibilities
    elseif vc==3
        V = inv(I)/n
    else
        V= Jinv*I*Jinv/n # sandwich form is preferred
    end
    return thetahat, objvalue, V, converged
end

function mle()
    println("mle(), with no arguments, runs mleresults(), which show usage")
    mleresults()
end
