using Statistics, LinearAlgebra, SolveDSGE

# this block reads and processes the file, leave it be
process_model("./CK.txt")
const dsge = retrieve_processed_model("./CK_processed.txt")

# solve model and simulate data
function dgp(θ, dsge, reps, rndseed=1234)
    p, ss = ParamsAndSS(θ)
    dsge = assign_parameters(dsge, p)
    scheme = PerturbationScheme(ss, 1.0, "third")
    solution = solve_model(dsge, scheme)
    burnin = 200
    nobs = 160
    data = simulate(solution, ss[1:3], reps*(burnin+nobs); rndseed = rndseed)
    # the next returns reps data sets, in an array of arrays
    data = [data[4:8, (nobs+burnin)*i-(nobs+burnin)+1+burnin:i*(nobs+burnin)]' for i = 1:reps]
end

function bad_data(data)
    any(isnan.(data)) || any(isinf.(data)) || any(std(data,dims=1) .==0.0) || any(data .< 0.0)
end

# this gives a vector of vectors, each a statistic drawn at the parameter value
function auxstat(θ, reps)
    auxstat.(dgp(θ, dsge, reps, rand(1:Int64(1e12))))
end

# These are the candidate auxiliary statistics for ABC estimation of
# the simple DSGE model of Creel and Kristensen (2013)
function auxstat(data)
    # check for nan, inf, no variation, or negative, all are reasons to reject
    if bad_data(data)
        return zeros(39)
    else    
        # known parameters
        α = 0.33
        δ = 0.025
        # recover capital (see notes)
        hours = data[:,3]
        intrate = data[:,4]
        wages = data[:,5]
        capital = α /(1.0-α )*hours.*wages./intrate
        # treat all variables
        logdata = log.([data capital])[2:end,:]
        # logs
        logoutput = logdata[:,1];   # output
        logcons = logdata[:,2];     # consumption
        loghours = logdata[:,3];    # hours
        logintrate = logdata[:,4];  # intrate
        logwages = logdata[:,5];    # wages
        logcapital = logdata[:,6]
        # rho1, sig1
        e = logoutput-α*logcapital-(1.0-α)*loghours 
        y = e[2:end]
        x = e[1:end-1]
        rho1 = cor(x,y)
        u = y-x*rho1
        sig1 = sqrt(u'*u/size(u,1))
        Z = vcat(rho1, sig1)
        # gam, rho2, sig2 (1/MRS=wage)
        x = [ones(size(logcons,1)) logcons]
        b = x\logwages
        e = logwages-x*b
        y = e[2:end]
        x = e[1:end-1]
        rho2 = cor(y,x)
        u = y-x*rho2
        sig2 = sqrt(u'*u/size(u,1))
        Z = vcat(Z, b, rho2, sig2)
        # standard devs. and correlations
        m = mean(logdata, dims=1)
        s = std(logdata, dims=1)
        d = (logdata .- m) ./s # keep means and std. devs., the VAR uses standardized and normalized   
        # AR(1)
        maxlag = 1
        y = d[2:end,:]
        x = d[1:end-1,:]
        n = size(y,1)
        rhos = zeros(6)
        es = zeros(n,6)
        for i = 1:6
            rho = x[:,i]\y[:,i]
            rhos[i] = rho
            es[:,i] = y[:,i]-rho.*x[:,i]
        end        
        varv = vech(cov(es)) # AR(1) error covariance elements 
        Z = vcat(Z, m[:], s[:], varv)
    end
    Z
end

function vech(x)
    k = size(x,1)
    a = zeros(Int((k^2-k)/2 + k))
    m = 1
    for i = 1:k
        for j = 1:i
            a[m] = x[i,j]
            m += 1
        end
    end
    a
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

function PriorSupport()
    lb = [0.95, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0/24.0]
    ub = [0.995, 5.0, 0.995, 0.1, 0.995, 0.1, 9.0/24.0]
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


