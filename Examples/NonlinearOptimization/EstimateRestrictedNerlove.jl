# Nerlove model by OLS (fminunc)
# and fmincon, imposing that factor shares are in (0,1) and sum to 1
using DelimitedFiles, Econometrics
# prepare the data
data = readdlm("nerlove.data")
data = log.(data[:,2:end])
n = size(data,1)
y = data[:,1]
x = [ones(n,1) data[:,2:5]]

# bounds and restriction
lb = [-1e6, -1e6, 0., 0., 0.0]
ub = [1e6, 1e6, 1., 1., 1.]
R = [0. 0. 1. 1. 1.]
r = 1.0

# define the objective function and start value
obj = theta -> (y-x*theta)'*(y-x*theta)
startval = (ub+lb)/2.0

# OLS
thetahat, objvalue = fminunc(obj, startval) 
println("the OLS estimates: obj. value: ", round(objvalue,digits=5))
prettyprint(thetahat)

# restricted LS
thetahat, objvalue_r, flag = fmincon(obj, startval, R, r, lb, ub) # both lower and upper bounds
println("the restricted LS estimates: obj. value: ", round(objvalue_r,digits=5))
prettyprint(thetahat)

println("Exercise: you should construct a qF test using the unrestricted and")
println("restricted objective function values")
