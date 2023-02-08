# plots the contours of the OLS objective function,
# with the true and estimated coefficients
# generate data
using Plots, Econometrics
n = 20
x = [ones(n) randn(n)]
beta = randn(2)
y = x*beta+randn(n)
# ols estimate and objective function
betahat = x\y
s = (a,b)-> (1/n)*sum((y - x*[a,b]).^2)
# plot the contours and the points
closeall()
b1 = range(-4.0,stop=4.0,length=100)
b2 = range(-4.0,stop=4.0,length=100)
p = contour(b1,b2,(b1,b2)->s(b1,b2),fill=true, c=:viridis)
scatter!([beta[1]],[beta[2]], markersize=10, label="true")
scatter!([betahat[1]],[betahat[2]], markersize=10, label="estimated")
# do OLS by simulated annealing, and plot the path
lb = -4*ones(2)
ub = -lb
thetastart = fill(-3.0, 2)
obj = theta-> (1/n)*sum((y - x*theta).^2)
results = samin(obj, thetastart, lb, ub, verbosity=2)
path = results[4]
plot!(path[:,4], path[:,5])
