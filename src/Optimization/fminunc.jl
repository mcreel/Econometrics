"""
    xopt, fopt, converged = fminunc(obj, startval)

Minimize the function obj, starting at startval. This is a convenience function, provided for 
students. I recommend learning to use Optim.jl rather than using this.

fminunc() with no arguments will run an example, execute edit(fminunc,()) to see the code.
"""
function fminunc(obj, x; xtol = 1e-08, gtol=1e-6, ftol=1e-12)
    results = try
        Optim.optimize(obj, x, LBFGS(), 
                            Optim.Options(
                            g_tol = gtol,
                            x_tol = xtol,
                            f_tol=ftol); autodiff=:forward)
    catch
        Optim.optimize(obj, x, LBFGS(), 
                            Optim.Options(
                            g_tol = gtol,
                            x_tol = xtol,
                            f_tol = ftol))
    end    
    return results.minimizer, results.minimum, Optim.converged(results)
end

function fminunc(silent=false)
    if !silent
        println("with no arguments, fminunc() runs a simple example")
        println("type edit(fminunc, ()) to see the example code")
    end
    # return of objective should be real valued, thus the [1] to pull value out of 1-dim array
    obj = x -> (x.^2)[1]
    # obj = x -> x'x # this would also work, as returns a Float64
    x = [2.0] # argument to objective function should be a vector, thus the brackets
    results = fminunc(obj, x)
    return results[3] # convergence check, for testing
end


