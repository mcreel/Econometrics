function ridge(y, x, k)
    # standardize regressors (except the constant)
    s = std(x,dims=1)
    A = eye(size(x,2))
    for i = 1:size(x,2)
        if s[i]==0
            s[i] = 1.0
            A[i,i] = 0.0 # don't shrink constant
        end
    end    
    x = x ./ s
    βhat = (inv(x'x + k*A)*x'y) ./ s'
end


