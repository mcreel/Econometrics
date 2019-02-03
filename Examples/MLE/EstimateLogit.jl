# Example of MLE estimation. The data is really drawn from the
# Logit DGP, so the model is well specified, and the MLE has
# the properties discussed in class. E.g., if you make n very
# large you should see that the estimator is very close to the
# true value of theta used to generate data

using Econometrics
include("LogitDGP.jl")
n = 30 # sample size
theta = [0, 0.5] # true theta for generating data
(y, x) = LogitDGP(n, theta) # generate the data

# now define things for estimation
model = theta -> logit(theta, y, x)
theta = zeros(size(x,2)) # start values for estimation

# Perform the estimation - Make sure that you examine
# the MLE estimation programs so that you see how this works
thetahat, objvalue, V, converged = mleresults(model, theta, "estimate logit model");

