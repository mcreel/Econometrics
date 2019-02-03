# this illustrates effect and detection of influential observations
using Plots, LinearAlgebra
n = 20
x = (1:n-1)/(n-1)
x = [x; 3] # the last observation is an outlying value of x

x = [ones(n,1) x]
P = x*inv(x'*x)*x'

beta = [10, -1]
e = 2*randn(n,1)
y = x*beta + e


# The fit
yhat = P*y

# calculate leverage and influence
leverage = diag(P)
e = y - yhat
influence = (leverage ./ (1.0 .-leverage)) .* e

xlabel = "X"
x = x[:,2]
scatter(x, y, label = "Data points")
plot!(x, yhat, label = "fitted")
plot!(x, leverage, label = "Leverage")
plot!(x, influence, label = "Influence")
gui()
#savefig("InfluentialObservation.svg")
