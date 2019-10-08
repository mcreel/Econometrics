# this is the log-likelihood for the linear model
# y = x*b+e, with e~N(0,s^2)
# the parameter theta = [b' s]
# this file is used for demonstration of MLE, to compare to OLS
function normal(theta, y, x)
    b = theta[1:end-1][1]
    s = theta[end][1]
    e = (y - x*b)./s
    logdensity = -log.(sqrt.(2.0*pi)) .- 0.5*log(s.^2) .- 0.5*e.*e
end
