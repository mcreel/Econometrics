% moments for OLS estimation
function [m, ms] = moments(theta, y, x)
    e = y - x*theta;
    n = size(y,1);
    k = size(theta,1);
    m = x'*e/n;
    if nargout > 1
        ms = x.*repmat(e,1,k);
    end    
end    
