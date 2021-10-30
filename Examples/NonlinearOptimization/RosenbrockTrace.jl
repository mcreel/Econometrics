using Optim, Plots
# method for contour
rosenbrock(x1,x2) = (1.0 - x1)^2 + 100.0 * (x2 - x1^2)^2
# method for optim
rosenbrock(θ) = (1.0 - θ[1])^2 + 100.0 * (θ[2] - θ[1]^2)^2
θ = zeros(2)
b1 = range(0.0,stop=1.5,length=100) 
b2 = range(0.0,stop=1.5,length=100)
contour(b1,b2,(b1,b2)->rosenbrock(b1,b2),fill=true, c=:viridis)
θs = zeros(100,2)
for i = 1:100
    #global θ = optimize(rosenbrock, θ, GradientDescent(), Optim.Options(iterations=1)).minimizer
    global θ = optimize(rosenbrock, θ, Newton(), Optim.Options(iterations=1)).minimizer
    θs[i,:] = θ'
    display(scatter!(θs[1:i,1], θs[1:i,2], legend=false))
    sleep(0.1)
end
θs


