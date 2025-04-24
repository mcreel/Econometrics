## Nerlove model by OLS (fminunc)
## and fmincon, imposing that factor shares are in (0,1) and sum to 1
using DelimitedFiles, Econometrics, Distributions
# prepare the data
cd(@__DIR__)
data = readdlm("../Data/nerlove.data")
data = log.(data[:,2:end])
n = size(data,1)
y = data[:,1]
x = [ones(n,1) data[:,2:5]]

# bounds and restriction
lb = [-1e6, -1e6, 0., 0., 0.] # factor shares non-negative
ub = [1e6, 1e6, 1., 1., 1.]   # factor shares not greater than 1
R = [0. 0. 1. 1. 1.] # homogeneity degree 1: factor shares sum to one
r = [1.0]

## define the objective function and start value
obj = theta -> (y-x*theta)'*(y-x*theta)
startval = (ub+lb)/2.0

## OLS
thetahat, objvalue = fminunc(obj, startval) 
println("the OLS estimates: obj. value: ", round(objvalue,digits=5))
prettyprint(thetahat)
# compare to analytic solution
x\y

## restricted LS: HOD1
thetahat, objvalue_r, flag = fmincon(obj, startval, R, r)
println("the restricted LS estimates: obj. value: ", round(objvalue_r,digits=5))
prettyprint(thetahat)

## restricted LS: HOD1 and shares in [0,1]
thetahat, objvalue_r, flag = fmincon(obj, startval, R, r, lb, ub)
println("the restricted LS estimates: obj. value: ", round(objvalue_r,digits=5))
prettyprint(thetahat)

println("Exercise: you should construct a qF test using the unrestricted and")
println("restricted objective function values")



