# compute R²
function rsq(Y, X)
    β = X\Y
    errors = Y - X*β
    R² = 1.0 - sum(errors.^2.0)/sum((Y .- mean(Y)).^2.0)
end
