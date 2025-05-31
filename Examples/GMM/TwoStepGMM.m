%%
NerloveData = importdata('NerloveData.m');
n = size(NerloveData,1);
y = log(NerloveData(:,2));
x = [ones(n,1) log(NerloveData(:,3:6))];
k = size(x,2);


%% this block shows that start values can be important, even for this simple problem
%% in the past, Matlab would not get the correct answer here, though Octave does. I'm
%% not sure of the current behavior of Matlab
%% Note April 2024: now Octave has copied Matlab's incorrect result! Bug for bug compatibility!
W = eye(k);
thetastart = zeros(k,1);
options = optimset('TolX',1e-15,'TolFun',1e-12);
[thetahat, obj] = fminunc(@(theta) moments(theta, y, x)'*W*moments(theta, y, x), thetastart);
fprintf('GMM results, identity weight matrix\n');
thetahat
obj


%% this finds the correct solution, using OLS start values (identical to GMM)
W = eye(k);
thetastart = x\y; % this give OLS coefficients
[thetahat, obj] = fminunc(@(theta) moments(theta, y, x)'*W*moments(theta, y, x), thetastart);
fprintf('GMM results, identity weight matrix\n');
thetahat
obj


%% this shows how to compute efficient weight matrix, and that it does
% not matter in this case (just identified)
[m ms] = moments(thetahat, y, x);
% compute efficient weight matrix
omega = ms'*ms/n;
fprintf('estimated covariance of moments\n');
omega
W = inv(omega);
thetastart = x\y; % this give OLS coefficients
[thetahat, obj] = fminunc(@(theta) moments(theta, y, x)'*W*moments(theta, y, x), thetastart);
fprintf('GMM results, optimal weight matrix\n');
thetahat
obj


%% how to compute covariance of estimates
% we're using the fact that moment contribs are independent here
% to replicate with GRETL, set cov estimator to HC0
D = -x'*x/n;
Vhat = (1/n)*inv(D*inv(omega)*D');
fprintf('estimated st. devs. of thetahat\n');
sqrt(diag(Vhat))


%% suppose we don't assume that moment contribs are independent
% then we should use Newey-West
lags = floor(size(x,1)^0.25); % choose lag length
omega = NeweyWest(ms, lags);
D = -x'*x/n;
Vhat = (1/n)*inv(D*inv(omega)*D');
fprintf('estimated st. devs. of thetahat, NW cov. estimator\n');
sqrt(diag(Vhat))
