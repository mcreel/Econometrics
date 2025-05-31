function omegahat = NeweyWest(Z,nlags)

% Returns the Newey-West estimator of the asymptotic variance matrix
% INPUTS: Z, a nxk matrix with rows the vector zt'
%         nlags, the number of lags
%
% OUTPUTS: omegahat, the Newey-West estimator of the covariance matrix


[n,k] = size(Z);

% de-mean the variables
Z = Z - repmat(mean(Z),n,1);

omegahat = Z'*Z/n; % sample variance
if nlags > 0
   % sample autocovariances
   for i = 1:nlags
      Zlag = Z(1:n-i,:);
      ZZ = Z(i+1:n,:);
      gamma = (ZZ'*Zlag)/n;
      weight = 1 - (i/(nlags+1));
      omegahat = omegahat + weight*(gamma + gamma');
   end
end



   
   











