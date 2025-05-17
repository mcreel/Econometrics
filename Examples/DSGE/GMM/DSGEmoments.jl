# computes the errors as outlined in the notes. Use these to generate
# moment conditions for estimation
function DSGEmoments(θ, data)
        # break out variables
        y = data[:,1]
        c = data[:,2]
        n = data[:,3]
        r = data[:,4]
        w = data[:,5]
        # break out params    
        α = 0.33
        δ = 0.025
        β, γ, ρ_z, σ_z, ρ_η, σ_η, nss = θ
        # recover ψ
        c1 = ((1.0/β + δ - 1.0)/α)^(1.0/(1.0-α))
        # compute ψ
        kss = nss/c1
        iss = δ*kss
        yss = kss^α * nss^(1-α)
        css = yss - iss
        ψ =  (css^(-γ)) * (1-α) * (kss^α) * (nss^(-α))
        # use MPL-MRS to get η, the preference shock: see Eqn. 16.17 in notes
        lnη = log.(w) .- log.(ψ) .- γ*log.(c)
        u = (lnη[2:end]-lnη[1:end-1]*ρ_η )/σ_η  # this should be N(0,1) 
        e1 = lnη[1:end-1].*u
        e2 = u.^2.0 .- 1.0 
        pref_shock = copy(u)
        # now the Euler eqn, see Eqn. 16.18 in notes
        euler = 100.0 *  ((1.0 .+ r .- δ).*β.*(c.^(-γ)) .- lag(c,1).^(-γ))[2:end]
        # recover K from MPK/MPL
        k = lag((α/(1.0-α))*n.*w./r,1)
        # get z, the production shock, from the production function
        z = log.(y) - α*log.(k) - (1.0-α)*log.(n)
        u = (z[2:end]-z[1:end-1]*ρ_z )/σ_z  
        e4 = z[1:end-1].*u
        e5 = u.^2.0 .- 1.0
        tech_shock = copy(u)
        # make moment conditions
        data, junk, junk = stnorm(log.(data)[1:end-1,:])
        errors = [100.0*e1 e2 100.0*e4 e5 euler data.*tech_shock data.*pref_shock data.*euler ]
        return errors
end

