# compute ols coefficients, fitted values, and errors
function lsfit(Y, X)
    beta = X\Y
    fit = X*beta
    errors = Y - fit
    return beta, fit, errors
end
