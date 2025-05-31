# uniform random walk, with bounds check
function proposal1(current, tuning)
    lb, ub = PriorSupport()
    trial = copy(current)
    if rand() > 0.0
        i = rand(1:size(current,1))
        trial[i] = current[i] + tuning[i].*randn()
    else
        trial = lb + (ub - lb).*rand(size(lb))
    end    
    return trial
end

function logL(θ, data)
    InSupport(θ) || return -Inf
    gs = DSGEmoments(θ, data)
    n, g = size(gs)
    Σ = cov(gs)
    isposdef(Σ) || return -Inf
    Σinv = inv(Σ)
    ghat = √n*mean(gs,dims=1)[:]
    lnL = - 0.5*dot(ghat, Σinv, ghat)
end


function logL_ABC(theta, data)
    n = size(data,1)
    gs = DSGEmoments(theta, data)
    Σ = cov(gs)
    cholΣ = chol(Σ)
    Σinv = inv(Σ)
    #Σinv = 1.0  # using the full asymptotic likelihood contracts too much
                # the true parameters can be outside the posterior's support
                # due to the small sample bias. Using a simple sum of squares
                # criterion doesn't have the problem
    ghat = mean(gs,1)
    #lnL = - 0.5*(ghat*Σinv*ghat')[1,1]
    lnL = log(sqrt(det(Σinv))) - 0.5*n*(ghat*Σinv*ghat')[1,1]
    ghat += randn(size(ghat))*cholΣ/sqrt(n) # used for nonparametric fitting
    return lnL, ghat


end


