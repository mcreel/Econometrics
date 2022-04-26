# Example of MLE estimation. The data is really drawn from the
# Logit DGP, so the model is well specified, and the MLE has
# the properties discussed in class. E.g., if you make n very
# large you should see that the estimator is very close to the
# true value of theta used to generate data

using Econometrics
include("LogitDGP.jl")
n = 30 # sample size
θ  = [0, 0.5] # true theta for generating data
(y, x) = LogitDGP(n, θ) # generate the data

# now define things for estimation
model = θ -> logit(θ, y, x)
θstart = zeros(size(x,2)) # start values for estimation

# Perform the estimation - Make sure that you examine
# the MLE estimation programs so that you see how this works
θhat, objvalue, V, converged = mleresults(model, θstart, "estimate logit model");

