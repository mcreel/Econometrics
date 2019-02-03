% MSM estimation for a sample from N(theta,1)
1;

function data = dgp(theta, n);
    data = theta + randn(n, 1);
end


function ms = msm_moments(theta, data)
    ms = zeros(size(data));
    n = 30;
    S = 10;
    for s = 1:S
        ms = ms + dgp(theta, n);
    end
    ms = data - ms/S;
end

% GMM
data = dgp(3, 30);
randn("seed",1234); % need to fix seed for simulated moments
W = 1;
thetahat = gmm_estimate(3, data, W, "msm_moments","");
thetahat
ms = msm_moments(thetahat, data);
W = inv(cov(ms));
gmm_results(thetahat, data, W, "msm_moments","");

