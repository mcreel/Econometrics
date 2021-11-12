# computes GARCH(1,1) log likelihood
using Statistics
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
