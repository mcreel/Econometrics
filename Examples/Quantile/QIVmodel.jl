using DelimitedFiles, Statistics
#=
 variables in card.data are
 log_wage
 educ
 exper
 black
 smsa
 south
 nearc4
 age
=# 
function getdata()
    data = readdlm("card.data")
    LNW = data[:,1]
    EDUC = data[:,2]
    EDUC = EDUC .- mean(EDUC)
    EXPER = data[:,3]
    EXPER = EXPER .- mean(EXPER)
    BLACK = data[:,4]
    SMSA = data[:,5]
    SOUTH = data[:,6]
    NEARC4 = data[:,7]
    AGE = data[:,8]
    EXPSQ = (EXPER.^2.0)
    constant = ones(size(LNW))
    Z = [constant NEARC4 AGE AGE.^2 BLACK SMSA SOUTH]
    X = [constant EDUC EXPER EXPSQ BLACK SMSA SOUTH]
    return LNW, X, Z
end

# the moments
function moments(beta::Array{Float64,1}, y::Array{Float64,1}, x::Array{Float64,2},
    z::Array{Float64,2}, tau::Float64)
    m = mean(z.*(tau .- (y .<= x*beta)),dims=1)
end

function likelihood(θ, y, x, z, τ, Σinv)
    m = moments(θ, y, x, z, τ)
    n = size(x,1)
    lnL = log(sqrt(det(Σinv))) + (-0.5*n*m*Σinv*m')[1]
end

function prior(theta)
    lb = [4.0, 0.0, 0.0, -0.1, -0.5, -0.5,  -0.5]
    ub = [9.0, 0.5, 0.2, 0.0, 0.5, 0.5,  0.5]
    a = 0.0
    if(all((theta .>= lb) .& (theta .<= ub)))
        a = 1.0
    end    
    return a
end

function proposal(current, tuning)
    i = rand(1:size(current,1))
    trial = copy(current)
    #trial[i] = trial[i] + tuning[i]*randn()
    trial = trial + tuning.*randn(size(trial))
    return trial
end

function proposal2(current, cholV)
    trial = copy(current)
    trial = trial + cholV'*randn(size(trial))
    return trial
end


