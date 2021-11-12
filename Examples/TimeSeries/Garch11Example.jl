# estimate GARCH(1,1) with AR(1) in mean
using Statistics, DelimitedFiles, Econometrics

# the likelihood function
function garch11(θ, y)
    # dissect the parameter vector
    μ, ρ, ω, α, β = θ
    ϵ = y[2:end] .- μ .- ρ*y[1:end-1]
    n = size(ϵ,1)
    h = zeros(n)
    # initialize variance; either of these next two are reasonable choices
    h[1] = var(y[1:10])
    #h[1] = var(y)
    for t = 2:n
        h[t] = ω + α*(ϵ[t-1])^2. + β*h[t-1]
    end
    logL = -log(sqrt(2.0*pi)) .- 0.5*log.(h) .- 0.5*(ϵ.^2.)./h
end

#function main()
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
#return
#end
#main()
