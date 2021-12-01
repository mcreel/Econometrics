# estimate GARCH(1,1) with AR(1) in mean
using Statistics, DelimitedFiles, Econometrics

# the likelihood function
function garch11(θ, y)
    # dissect the parameter vector
    μ, ρ, ω, α, β = θ
    ylag = y[1:end-1]
    y = y[2:end]
    ϵ = y .- μ .- ρ*ylag
    n = size(ϵ,1)
    h = zeros(n)
    # initialize variance; either of these next two are reasonable choices
    #h[1] = var(y[1:10])
    h[1] = var(y)
    for t = 2:n
        h[t] = ω + α*(ϵ[t-1])^2. + β*h[t-1]
    end
    logL = -log(sqrt(2.0*pi)) .- 0.5*log.(h) .- 0.5*(ϵ.^2.)./h
end

function main()
# weekly close price of NSYE, data provided with GRETL
data = readdlm("nysewk.txt")
# compute weekly percentage growth
y = 100.0 * log.(data[2:end] ./ data[1:end-1])
# Constrained maximization of logL
# note that the objective has a minus sign in front, as fmincon
# minimizes, but we want to maximize the logL
# get constrained estimates to use as input for MLE
θstart = [mean(y); 0.0; var(y); 0.1; 0.1]
obj = θ -> -mean(garch11(θ, y))
θhat, logL, junk  = fmincon(obj, θstart, [], [], [-Inf, -1.0, 0.0, 0.0, 0.0], [Inf, 1.0, Inf, 1.0, 1.0])
# use this to report results. Only works because constaints not bining.
model = θ -> garch11(θ, y)
mleresults(model, θhat, "Garch(1,1) results", ["μ", "ρ","ω","α","β"])
nothing
end
main()
