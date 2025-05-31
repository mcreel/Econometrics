using SolveDSGE

nfeatures = 5
nparams = 7

function PriorSupport()
    Float32.([0.95, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0/24.0]),      # lb
    Float32.([0.995, 5.0, 0.995, 0.1, 0.995, 0.1, 9.0/24.0])  # ub
end

function StNormParams()
    lb, ub = PriorSupport()
    (ub+lb)/2, sqrt.(((ub-lb).^2)/12)
end    

function InSupport(θ)
    lb, ub = PriorSupport()
    all(θ .< ub) && all(θ .>= lb)
end    

function TrueParameters()
 [
 0.99,  # β
 2.0,   # γ     
 0.9,   # ρ₁  
 0.02,  # σ₁   
 0.7,   # ρ₂  
 0.01,  # σ₂   
 8.0/24.0]  # nss
end    



function PriorDraw(S)
    lb, ub = PriorSupport()
    θ = zeros(Float32, 7, S)
    Threads.@threads for s ∈ axes(θ, 2)
        θ[:, s] = (ub - lb) .* rand(Float32,7) + lb
    end
    θ
end

function ParamsAndSS(params)
    α = 0.33
    δ = 0.025
    β, γ, ρ₁, σ₁, ρ₂, σ₂, nss = params
    c1 = ((1/β  + δ - 1)/α)^(1/(1-α))
    kss = nss/c1
    iss = δ*kss
    yss = kss^α * nss^(1-α)
    css = yss - iss;
    MUCss = css^(-γ)
    rss = α * kss^(α-1) * nss^(1-α)
    wss = (1-α)* (kss)^α * nss^(-α)
    MULss = wss*MUCss
    ψ =  (css^(-γ)) * (1-α) * (kss^α) * (nss^(-α))
    p = [β, γ, ρ₁ , σ₁, ρ₂, σ₂, ψ]
    ss = [0.0, 0.0, kss, yss, css, nss, rss, wss, MUCss, MULss]
    return p, ss
end   


# S simulations, each at a different θ, for training the net
# the training uses a second order perturbation, for speed
# θs are returned standardized and normalized, for training
function MakeData(S, model)
    θs = zeros(Float32, 7, S)
    nobs = 160
    X = zeros(Float32, 160, S, 5)
    # solve models and simulate data
    Threads.@threads for s = 1:S
        ok = false
        while !ok
            θ = PriorDraw(1)
            p, ss = ParamsAndSS(θ)
            model2 = assign_parameters(model, p)
            scheme = PerturbationScheme(ss, 1.0, "second")
            solution = solve_model(model2, scheme)
            burnin = 100
            rndseed = rand(1:Int64(1e10))
            data = simulate(solution, copy(ss[1:3]), burnin+nobs; seed = rndseed)[4:8, burnin+1:end]'
            ok = all(data .< 10.0) && all(data .> 0.0)
            if ok
                data .-= [0.84, 0.69, 0.33, 0.05, 1.72]'
                data ./= [0.51, 0.44, 0.36, 0.018, 0.34]'
                X[:, s, :] = Float32.(data)
                θs[:,s] = Float32.(θ)
            end
        end    
    end    
    tabular2conv(permutedims(Float32.(X), (3, 2, 1))), Float32.(TransformParameters(θs))
end

function TransformParameters(θ)
    m, s = StNormParams()
    (θ .- m) ./ s
end    

function UntransformParameters(θ)
    m, s = StNormParams()
    s.*θ .+ m
end    


# Want S simulations, all at the same θ
# This does one long simulation, and separates samples using burnin draws
# between them
function MakeData(θ, S, model)
    InSupport(θ) || throw(ArgumentError("θ is not in support"))
    nobs = 160
    X = zeros(Float32, nobs, S, 5)
    # solve models and simulate data
    p, ss = ParamsAndSS(θ)
    model = assign_parameters(model, p)
    scheme = PerturbationScheme(ss, 1.0, "second")
    solution = solve_model(model, scheme)
    burnin = 500
    rndseed = rand(1:Int64(1e10))
    data = simulate(solution, copy(ss[1:3]), S*(burnin+nobs); seed = rndseed)[4:8,:]'
    data = min.(data, 10.0)
    data = max.(data, 0.0)
    data .-= [0.84, 0.69, 0.33, 0.05, 1.72]'
    data ./= [0.51, 0.44, 0.36, 0.018, 0.34]'
     # split up the samples
    Threads.@threads for s ∈ axes(X, 2)
        X[:, s, :] = data[(nobs+burnin)*s-(nobs+burnin)+1+burnin:s*(nobs+burnin),:]
    end
    tabular2conv(permutedims(Float32.(X), (3, 2, 1)))
end



