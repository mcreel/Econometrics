# uniform random walk, with bounds check
function proposal1(current, tuning)
    lb_param_ub = [
        0.95    0.99   	0.999; 	    # beta
        0.0	    2      	5;		    # gam
        0    	0.9   	0.999;	    # rho1
        0.0001       0.02 	0.1;    # sigma1
        0    	0.7     0.999;      # rho2
        0.0001	    0.01  	0.1;    # sigma2
        6/24    8/24	9/24	    # nss
    ]
    lb = lb_param_ub[:,1]
    ub = lb_param_ub[:,3]
    trial = copy(current)
    if rand() > 0.0
        i = rand(1:size(current,1))
        trial[i] = current[i] + tuning[i].*randn()
    else
        trial = lb + (ub - lb).*rand(size(lb))
    end    
    return trial
end

function proposal2(current, cholV)
    trial = copy(current)
    lb_param_ub = [
        0.95    0.99   	0.999; 	    # beta
        0.0	    2      	5;		    # gam
        0    	0.9   	0.999;	    # rho1
        0.0001       0.02 	0.1;    # sigma1
        0    	0.7     0.999;      # rho2
        0.0001	    0.01  	0.1;    # sigma2
        6/24    8/24	9/24	    # nss
    ]
    lb = lb_param_ub[:,1]
    ub = lb_param_ub[:,3]
    if rand() > 0.0
        trial += cholV'*randn(size(trial))
    else
        trial = lb + (ub - lb).*rand(size(lb))
    end


end

function prior(theta)
    lb_param_ub = [
        0.95    0.99   	0.999; 	# beta
        0.0	    2      	5;		    # gam
        0    	0.9   	0.999;	    # rho1
        0.0001       0.02 	0.1;    # sigma1
        0    	0.7     0.999;     # rho2
        0.0001	    0.01  	0.1;    # sigma2
        6/24    8/24	9/24	    # nss
    ]
    lb = lb_param_ub[:,1]
    ub = lb_param_ub[:,3]
    a = 0.0
    if(all((theta .>= lb) .& (theta .<= ub)))
        a = 1.0
    end
    return a
end

function logL(theta, data)
    n = size(data,1)
    gs = DSGEmoments(theta, data)
    sigma = cov(gs)
    siginv = inv(sigma)
    ghat = mean(gs,dims=1)
    #lnL = - 0.5*(ghat*siginv*ghat')[1,1]
    lnL = log(sqrt(det(siginv))) - 0.5*n*(ghat*siginv*ghat')[1,1]
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


