# computes the errors as outlined in the notes. Use these to generate
# moment conditions for estimation
function DSGEmoments(thetas, data)
        # break out variables
        y = data[:,1]
        c = data[:,2]
        n = data[:,3]
        r = data[:,4]
        w = data[:,5]
        # break out params    
        alpha = 0.33
        delta = 0.025
        beta = thetas[1]
        gam = thetas[2]
        rho_z = thetas[3]
        sig_z = thetas[4]
        rho_eta = thetas[5]
        sig_eta = thetas[6]
        nss = thetas[7]
        # recover psi
        c1 = ((1.0/beta + delta - 1.0)/alpha)^(1.0/(1.0-alpha))
        kss = nss/c1
        iss = delta*kss
        yss = kss^alpha * nss^(1-alpha)
        css = yss - iss
        psi =  (css^(-gam)) * (1-alpha) * (kss^alpha) * (nss^(-alpha))
        # use MPL-MRS to get eta, the preference shock
        eta = log.(w) .- gam*log.(c) .- log.(psi)
        X = lag(eta,1)
        u = (eta-X*rho_eta)/sig_eta
        e1 = X.*u
        e2 = u.^2.0 .- 1.0
        pref_shock = copy(u)
        # now the Euler eqn
        e3 = (1.0 .+ r .- delta).*beta.*(c.^(-gam)) .- lag(c,1).^(-gam) 
        e3 = 100.0*e3
        # recover K from MPK/MPL
        lagk = (alpha/(1.0-alpha))*lag(n.*w./r,1)
        # get z, the production shock, from the production function
        z = log.(y) - alpha*log.(lagk) - (1.0-alpha)*log.(n)
        X = lag(z,1)
        u = (z-X*rho_z)/sig_z
        e4 = X.*u
        e5 = u.^2.0 .- 1.0
        tech_shock = copy(u)
        # make moment conditions
        data = lag(data,1)
        errors = [e1 e2 e3 e4 e5 pref_shock.*data tech_shock.*data data.*e3]
        errors = errors[2:end,:] # need to drop 2, because lagk uses a lag, and we use lagged lagk
        #dstats(errors);
        return errors
end

