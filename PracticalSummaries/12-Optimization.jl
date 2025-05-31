# Practical Review: Numeric Optimization (Ch. 12)

##
# Let's start with a simple problem, where we already know the
# answers. We will apply numeric optimization methods and we 
# will be able to verify that they give use the correct results.
# To start, let's define some data, and do the basic OLS 
# computations
using Econometrics
n = 5
x = [ones(n) 1:n]
β = [10., -1.]
ϵ = randn(n)
y = x*β + ϵ
ols(y,x)

##
# Let's see how to get these results using numeric minimization
# and extremum estimation theory...
using Econometrics, Statistics
objᵢ = β -> (y-x*β).^2.  	# iᵗʰ obs contrib to sum of squares
obj = β -> mean(objᵢ(β))	# average objective
βhat, objval = fminunc(obj, zeros(2))

##
# fminunc is just a frontend to the Optim.jl package's LBFGS unconstrained
# minimizer routine. To call that directly, which exposes more options, 
# such as choosing the algorithm, do
using Optim
tol = 1e-08  # you want to use tight enough tolerances to ensure an accurate answer
results = Optim.optimize(obj, zeros(2), LBFGS(), 
                            Optim.Options(
                            g_tol = tol,
                            x_tol=tol,
                            f_tol=tol); autodiff=:forward)

##
# The output has the function value and convergence information. To see the
# value of the minimizer, do
βhat = results.minimizer

##
# Let's check that the Hessian is p.d.
using ForwardDiff, LinearAlgebra
H = ForwardDiff.hessian(obj, βhat)
eigen(H).values
# The eigenvalues are all positive, so the Hessian is p.d. at
# the solution. For the OLS case, the objective function is
# known to be globally convex, so multiple local optima are not
# a problem.

## Let's look at another problem, which has multiple
# minima, and see how to find the global minimizer
# Here's another objective function
f(x,y) = -(3*(1-x)^2 * exp(-(x^2) - (y+1)^2)
    - 6*(x/5 - x^3 - y^5) * exp(-x^2 - y^2)
    - 1/3 * exp(-(x+1)^2 - y^2))
# let's explore it:
using Plots
plotly()
x = range(-4, step=0.1, stop=4)
y = x
surface(x, y, (x,y)->f(x,y), color=:viridis)

##
contour(x,y,(x,y)->f(x,y), c=:viridis)
xlabel!("x")
ylabel!("y")

##
# This function has 3 local min, and 3 local max.
# The global minimizer is near (0,1.6), with a function
# value of -4.86. We will need some care to find the global
# minimizer. The domain is (-4,4)×(-4,4), which may help us.
# First, let's try gradient-based minimization, starting at
# a random point in the domain.

##
# Before we do that, the form of the function defined above
# does not take a single vector argument. For optimization, we
# need all parameters in a single vector. Julia's multiple
# dispatch allow us to do this easily, by defining another
# method for the function f:
function f(θ)
   x,y = θ
   f(x,y) # pretty neat, eh?
end

##
#= 
Another thing to check before doing optimization! Remember
that gradient based methods will work better when the data
are scaled to put the elements of the gradient on a common
order of magnitude. You can check this by examining the
gradient at some points. Run this block a few times to check.
=#
using ForwardDiff
lower = [-4., -4.]
upper = -lower
θstart = 8.0 .* rand(2) .-4.0 
ForwardDiff.gradient(f, θstart)
# in this case, things seem to be ok.

##
# Now we can try optimization. Running this block repeatedly (Alt-enter)
# will show that each of the 3 minimizers are found. Thus, if one defined
# LBFGS with enough start values, the global minimizer would be found.
# If you run the code enough times, you will see that the algorithm
# occasionally goes off far away from the 3 local minimizers. In that
# region, the gradients are zero, as we see in the plots, and the algorithm
# can't find better points.

# NOTE: in class, split the plot pane out to a new window, to do multiple runs

using Econometrics, Plots, Optim
gr() # using the GR backend
x = range(-4, step=0.1, stop=4)
y = x
Plots.contour(x,y,f,c=:viridis)
Plots.xlabel!("x")
Plots.ylabel!("y")

θstart = 8.0 .* rand(2) .-4.0 
tol = 1e-08
results = Optim.optimize(f, θstart, LBFGS(), 
                            Optim.Options(
                            g_tol = tol,
                            x_tol=tol,
                            f_tol=tol); autodiff=:forward)

θhat = results.minimizer
Plots.scatter!([θhat[1]],[θhat[2]],legend=false, markersize=10, c = :viridis)

##
# A more robust algorithm is to use the known domain to define
# box constraints. This will still bounce between the 3 local
# min, but it will not shoot off to far away points.
using Econometrics, Plots, Optim
gr() # using the GR backend
x = range(-4, step=0.1, stop=4)
y = x
Plots.contour(x,y,f, c=:viridis)
Plots.xlabel!("x")
Plots.ylabel!("y")

lower = [-4., -4.]
upper = -lower
θstart = 8.0 .* rand(2) .-4.0 
tol = 1e-08
results = Optim.optimize(f, lower, upper, θstart, Fminbox(LBFGS()), 
                            Optim.Options(
                            g_tol = tol,
                            x_tol=tol,
                            f_tol=tol); autodiff=:forward)
θhat = results.minimizer
Plots.scatter!([θhat[1]],[θhat[2]],legend=false, markersize=10)

## 
#
# Global optimization
#
# Another option is to use a global optimizer. An example is the
# simulated annealing algorithm. This one should find the global
# min every time, provided the cooling schedule is slow enough
using Econometrics, Plots
# try setting rt=0.05, it will fail sometimes, because the 
# search region is shrunk too quickly
# rt=0.9 has worked every time I've tried it
plotly() # using the GR backend
x = range(-4, step=0.1, stop=4)
y = x
Plots.contour(x,y,f, c=:viridis)

rt = 0.8  # This controls the cooling schedule. 
            # If close enough to 1, gives a good guarantee.
            # Closer to 0 is faster, but more likely to miss the 
            # global minimizer
lower = [-4., -4.]
upper = -lower
θstart = [2.0, -3.9] # a pretty bad starting point
trace = (samin(f, θstart, lower, upper, rt=rt)[4])[:,4:5]
# plot the path take to the solution (run several times to see examples)
for i = 1:size(trace,1)-1
    display(Plots.scatter!([trace[i,1]], [trace[i,2]], labels=:none, color=:red, markersize=3))
    sleep(0.05)
end
θhat = round.(trace[end,:], digits=5)
display(scatter!([trace[end,1]], [trace[end,2]], legend=true, label="θhat=$θhat", color=:yellow, markershape=:circle, markersize=5))
println("for reference, the global minimizer is (-0.00879, 1.58163)")
println("with the function value -4.8611920055")

#=
To summarize, a gradient based minimizer should be reliable
when the objective function is globally convex. When this is
not the case, one can use multiple start values for a gradient
based method, or one might prefer to use a global minimizer.
=#

