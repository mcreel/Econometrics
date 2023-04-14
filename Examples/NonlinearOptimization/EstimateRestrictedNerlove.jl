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
lb = [-1e6, -1e6, 0., 0., 0.0]
ub = [1e6, 1e6, 1., 1., 1.]
R = [0. 0. 1. 1. 1.]
r = [1.0]

## define the objective function and start value
obj = theta -> (y-x*theta)'*(y-x*theta)
startval = (ub+lb)/2.0

## OLS
thetahat, objvalue = fminunc(obj, startval) 
println("the OLS estimates: obj. value: ", round(objvalue,digits=5))
prettyprint(thetahat)

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
println("Because of inequality restrictions, the asyptotic distribution of the test is complicated.
Nevertheless, just using the critical value from the regular χ² distribution suggests that the
restrictions are not rejected at conventional significance levels.")


