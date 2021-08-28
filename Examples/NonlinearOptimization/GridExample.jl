using Econometrics, Plots
function GridExample(points, doprint=false)
    # plot the line
    x = range(-pi,stop=pi/2.0,length=1000)
    obj = theta-> 0.5*theta*sin(theta^2.0)
    y = obj.(x)
    plot(x, y, legend=false)
    # plot the grid points
    x = range(-pi,stop=pi/2.0,length=points)
    y = obj.(x)
    scatter!(x, y)
    # plot the best point found
    scatter!([x[argmin(y)]], [y[argmin(y)]], markersize=10)
    p = plot!(legend=false)
    if doprint savefig("gridsearch.png") end
    return p
end
