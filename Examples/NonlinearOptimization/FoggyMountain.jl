using Econometrics
include("FoggyMountainObj.jl")

function main()
## One BFGS run with poor starting values
θ = [8.0, -8.0]
θhat, obj_value, converged = fminunc(FoggyMountainObj, θ)
println()
println("The result with poor start values")
println("objective function value: ", round(obj_value,digits=5))
println("local minimizer: ", θhat)
println()

## Now try simulated annealing
lb = [-20.0,-20.0]
ub = -lb
xopt = samin(FoggyMountainObj, θ, lb, ub, verbosity=1)


## Now try BFGS with multiple start values
bestobj = 1000
θhat = zeros(2)
for i = 1:50
    θ = (ub-lb).*rand(2) + lb # use same bounds as SA
    θ, obj_value, converged = fminunc(FoggyMountainObj, θ)
    if obj_value < bestobj
        θhat = θ
        bestobj = obj_value
    end
end
println("the estimate using BFGS with multiple start values: ", θhat)
println("the best objective function value: ", bestobj)
end
nothing
main()
