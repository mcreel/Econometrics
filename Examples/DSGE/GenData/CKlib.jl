# this block reads and processes the file, leave it be
filename = "CK.txt"
path = joinpath(@__DIR__,filename)
process_model(path)
processed_filename = "CK_processed.txt"
processed_path =  joinpath(@__DIR__,processed_filename)
dsge = retrieve_processed_model(processed_path)

function CKdgp(θ, dsge, rndseed=1234)
    p, ss = ParamsAndSS(θ)
    dsge = assign_parameters(dsge, p)
    scheme = PerturbationScheme(ss, 1.0, "third")
    solution = solve_model(dsge, scheme)
    burnin = 1000
    nobs = 160
    data = simulate(solution, ss[1:3], burnin+nobs; rndseed = rndseed)
    data = data[4:8, burnin+1:end]'
end

dgp = θ -> CKdgp(θ, dsge, rand(1:Int64(1e12)))

function TrueParameters()
[0.33,  # α 
 0.99,  # β
 0.025, # δ 
 2.0,   # γ     
 0.9,   # ρ₁  
 0.02,  # σ₁   
 0.7,   # ρ₂  
 0.01,  # σ₂   
 8.0/24.0]  # nss
end    

function PriorSupport()
    lb = [0.33, 0.95, 0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0/24.0]
    ub = [0.33, 0.995, 0.025, 5.0, 0.995, 0.1, 0.995, 0.1, 9.0/24.0]
    lb,ub
end    

function PriorDraw()
    lb, ub = PriorSupport()
    θ = (ub-lb).*rand(size(lb,1)) + lb
end    

function InSupport(θ)
    lb,ub = PriorSupport()
    all(θ .>= lb) & all(θ .<= ub)
end

function Prior(θ)
    InSupport(θ) ? 1.0 : 0.0
end    


function ParamsAndSS(params)
    α, β, δ, γ, ρ₁, σ₁, ρ₂, σ₂, nss = params
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
    p = [α, β, δ, γ,ρ₁ , σ₁, ρ₂, σ₂, ψ]
    ss = [0.0, 0.0, kss, yss, css, nss, rss, wss, MUCss, MULss]
    return p, ss
end   


