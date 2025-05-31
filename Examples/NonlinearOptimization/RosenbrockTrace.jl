using Optim, Plots
gr()
# method for contour
rosenbrock(x1,x2) = (1.0 - x1)^2 + 100.0 * (x2 - x1^2)^2
# method for optim
rosenbrock(θ) = (1.0 - θ[1])^2 + 100.0 * (θ[2] - θ[1]^2)^2
θ = zeros(2)
b1 = range(0.0,stop=1.5,length=100) 
b2 = range(0.0,stop=1.5,length=100)
contour(b1,b2,(b1,b2)->rosenbrock(b1,b2),fill=true, c=:viridis)
iters = 50
θs = zeros(iters,2)
# plot the true minimizer
display(scatter!([1.0, 1.0], color=:yellow, markershape=:star5, markersize=20))
for i = 1:iters    
    #global θ = optimize(rosenbrock, θ, GradientDescent(), Optim.Options(iterations=1)).minimizer
    global θ = optimize(rosenbrock, θ, Newton(), Optim.Options(iterations=1)).minimizer
    θs[i,:] = θ'
    display(scatter!(θs[1:i,1], θs[1:i,2], legend=false))
    gui()
    sleep(0.05)
end

θs


