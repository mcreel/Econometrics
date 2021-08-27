# example of CUE GMM: draws from N(0,1)
using Econometrics, Statistics, LinearAlgebra
y = randn(100)
# 3 moment conditions
moments = theta -> [(y .-theta[1]) ((y.^2.0) .- theta[2]) ((y .- theta[1]).^3.0)]
m = theta -> mean(moments(theta), dims=1)[:]
weight = theta -> inv(NeweyWest(moments(theta)))
obj = theta -> dot(m(theta), weight(theta), m(theta))
theta = [0.0, 1.0]
results = fminunc(obj, theta)


