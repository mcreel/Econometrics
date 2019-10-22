using Distributions, Plots, Random
"""
simple univariate kernel density that uses rule-of-thumb
bandwidth and Gaussian kernel
* x is n-vector of data
* xeval is neval-vector of evaluation points
* output is neval-vector of estimated densities

execute ndens() for an example
"""
function npdens()
    println("npdens(), with no arguments, runs a simple example")
    println("execute edit(ndens,()) to see the code")
    Random.seed!(1) # set seed to enable testing
    bandwidth = 1.0
    x = rand(Chisq(5),1000)
    xeval, d = npdens(x)
    dt = pdf.(Ref(Chisq(5)), xeval)
    Plots.plot(xeval, d, label="kernel")
    Plots.plot!(xeval, dt, label="true")
end 

function npdens(x::Array{Float64,1}; points::Int=1000)
    xeval = Array(range(minimum(x), stop=maximum(x), length = points))
    dens = zero.(similar(xeval))
    n = size(x,1)
    bandwidth = 1.06*std(x)/(n^0.2) # rule-of-thumb: replace with CV?
    dist = Normal(0.0,1.0)
    for i = 1:points
        dens[i] = mean(pdf.(Ref(dist), (x .- xeval[i])./bandwidth))
    end
    dens ./= bandwidth
    return xeval, dens
end
