# estimate GARCH(1,1) with AR(1) in mean
using Statistics, DelimitedFiles, Econometrics
function main()
include("garch11.jl")
# weekly close price of NSYE, data provided with GRETL
data = readdlm("nysewk.txt")
# compute weekly percentage growth
y = 100.0 * log.(data[2:end] ./ data[1:end-1])
# Constrained maximization of logL
# note that the objective has a minus sign in front, as fmincon
# minimizes, but we want to maximize the logL
# get constrained estimates to use as input for MLE
thetastart = [mean(y); 0.0; var(y); 0.1; 0.1]
obj = theta -> -sum(garch11(theta, y))
thetahat, logL, junk  = fmincon(obj, thetastart, [], [], [-Inf, -1.0, 0.0, 0.0, 0.0], [Inf, 1.0, Inf, 1.0, 1.0])
# now do MLE to see formatted results.
# NOTE TO SELF: this won't work when the constraints are binding
# should add a method to mle for constrained problems.
model = theta -> garch11(theta, y)
thetahat, logL, junk, converged = mleresults(model, thetahat, "GARCH(1,1) example")
return
end
main()
