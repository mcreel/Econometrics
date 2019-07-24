# compute R²
function rsq(Y, X)
    β = X\Y
    errors = Y - X*β
    R² = 1.0 - sum(errors.^2.0)/(size(X,1)-size(X,2))
end
