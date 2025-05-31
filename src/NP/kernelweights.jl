using Distances
# kernelweights.jl: kernel weights for nonparametrics
function kernelweights(x, xeval, bandwidth, prewhiten=true, kernel="gaussian", neighbors=1)
    if prewhiten
        v = cov(x)
        if isposdef(v)
            P = inv((cholesky(v).U))
        else
            P = sqrt(inv(diag(v)))
        end    
        x *= P
        xeval *= P
    end    
    neval = size(xeval,1)
    n = size(x,1)
    weights = zeros(n,neval)
    # Gaussian product kernel, add others if desired
    if kernel=="gaussian"
        for i=1:neval
            z = (x.-xeval[i:i,:])/bandwidth
            weights[:,i] = exp.(-0.5*sum(z.*z,dims=2))
        end
    end

    if kernel=="knngaussian" 
        distances = pairwise(Euclidean(),x', xeval', dims=2) # get all distances
        @inbounds  for i = 1:neval
            di = distances[:,i]
            ind = sortperm(di) # indices of k nearest neighbors
            selected = vec(ind[1:neighbors,:])
            z = (x[selected,:].-xeval[[i],:])/bandwidth
            weights[selected,i] = exp.(-0.5*sum(z.*z,dims=2))
        end
    end    
    weights ./= sum(weights,dims=1)
    return weights
end
