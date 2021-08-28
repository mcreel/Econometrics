include(@__DIR__*"FoggyMountainObj.jl")

function main()
# One BFGS run with poor starting values
theta = [8.0, -8.0]
thetahat, obj_value, converged = fminunc(FoggyMountainObj, theta)
println()
println("The result with poor start values")
println("objective function value: ", round(obj_value,digits=5))
println("local minimizer: ", thetahat)
println()

# Now try simulated annealing
lb = [-20.0,-20.0]
ub = -lb
xopt = samin(FoggyMountainObj, theta, lb, ub, verbosity=1)

# Now try BFGS with multiple start values
bestobj = 1000
thetahat = zeros(2)
for i = 1:50
    theta = (ub-lb).*rand(2) + lb # use same bounds as SA
    theta, obj_value, converged = fminunc(FoggyMountainObj, theta)
    if obj_value < bestobj
        thetahat = theta
        bestobj = obj_value
    end
end
println("the estimate using BFGS with multiple start values: ", thetahat)
println("the best objective function value: ", bestobj)
end
main()
