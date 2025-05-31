# Estimates the basic Nerlove Cobb-Douglas model by MLE
# this is intended to show that MLE with normality is the
# ML estimator (for the b"s, not for sig^2)
using Econometrics, DelimitedFiles
cd(@__DIR__)
data = readdlm("../Data/nerlove.data")
data = data[:,2:end]
data = log.(data)
n = size(data,1)
y = data[:,1]
x = data[:,2:end]
x = [ones(n,1) x]
names = ["constant", "output", "labor", "fuel", "capital", "sig"]
theta = [zeros(size(x,2)); 1.0] # start values for estimation
model = theta -> normal(theta, y, x)
thetahat, objvalue, V, converged = mleresults(model, theta, "estimate Nerlove model by MLE", names)
nothing

