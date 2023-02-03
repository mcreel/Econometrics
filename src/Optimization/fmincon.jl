"""
    xopt, fopt, converged = fmincon(obj, startval, R, r, lb, ub; tol=tol, iterlim=0)

Minimize the function obj subject to linear restraints
Rb=r and box constraints defined by lower bounds lb 
and upper bounds ub, starting at startval. 

R and r can both be empty, by passing []
lb and/or ub can also be empty, by passing []

function and parameter tolerance can be set using keyword tol (default is 1e-10)
iteration limit can be set with keyword iterlim (default is no limit)

fmincon() with no arguments will run an example, execute edit(fmincon,()) to see the code.
fmincon() uses NLopt.jl to do the actual minimization.

"""
function fmincon(obj, startval, R=[], r=[], lb=[], ub=[]; tol = 1e-10, iterlim=0)
    # the objective is an anonymous function
    function objective_function(x::Vector{Float64}, grad::Vector{Float64})
        obj_func_value = obj(x)[1,1]
        return(obj_func_value)
    end
    # impose the linear restrictions
    function constraint_function(x::Vector, grad::Vector, R, r)
        result = R*x .- r
        return result[1,1]
    end
    opt = Opt(:LN_COBYLA, size(startval,1))
    min_objective!(opt, objective_function)
    # impose lower and/or upper bounds
    if lb != [] lower_bounds!(opt, lb) end
    if ub != [] upper_bounds!(opt, ub) end
    # impose linear restrictions, by looping over the rows
    if R != []
        for i = 1:size(R,1)
            equality_constraint!(opt, (theta, g) -> constraint_function(theta, g, R[i:i,:], r[i]), tol)
        end
    end    
    xtol_rel!(opt, tol)
    ftol_rel!(opt, tol)
    maxeval!(opt, iterlim)
    (objvalue, xopt, flag) = NLopt.optimize(opt, startval)
    return xopt, objvalue, flag
end

function fmincon(silent=false)
    if !silent
        println("with no arguments, fmincon() runs a simple example")
        println("type edit(fmincon, ()) to see the example code")
    end
    # return of objective should be real valued, thus the [1] to pull value out of 1-dim array
    obj = x -> x'x
    x = [2.0, 2.0]
    # sum of params should be 1
    R = [1.0 1.0]
    r = 1.0
    results = fmincon(obj, x, R, r)
    # check convergence
    results[3] == :FTOL_REACHED
end





