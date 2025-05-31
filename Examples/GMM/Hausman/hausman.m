# Hausman test comparing OLS and GIV estimators
# ref. Creel (2004) "Modified Hausman tests for inefficient estimators" Applied Economics
# this calculates the H2 test statistic, which uses the joint moment conditions and
# the overall efficient joint weight matrix

# parameters of DGP (see the function just below)
# you can vary these, but not that there are limits
# on the values of rho and good_inst that will lead
# to a positive definite covariance. You'll get an
# error message if you set them too high.

n = 200; # sample size
good_inst = 0.5; # correlation btw. regressor and instrument
rho = 0.0; # correlation btw. error and regressor
# when rho is not zero, the test should start to show power. When rho=0, the percentage
# of times that the p-value is less than a given significance level should be equal to
# the significance level. You can Monte Carlo this to verify the good behavior of the
# test statistic.

# simple dgp with a single regressor and instrument
# Note that test results are invariant to any scalar multiple
# of the variance of the error, so it's set to one.
#
# good inst is the correlation between the regressor and the instrument
# rho is the correlation between the regressor and the error
function [y, x, w] = dgp_linear(n, good_inst, rho)

	sig = [1,        rho,  good_inst;
		   rho,       1,    0        ;
		   good_inst, 0,    1       ];
		
	vars = randn(n,3) * chol(sig);
	e = vars(:,2);
	w = vars(:,3);
	x = vars(:,1);
	e = e .* (ones(n,1) +  x); % heteroscedastic: OLS is not efficient
	y = x + e;
endfunction


# moments for both OLS and IV
function m = joint(theta, data, momentargs)
	k = momentargs{1};
	theta1 = theta(1:k,:);
	theta2 = theta(k+1:2*k,:);
	m1 = qml(theta1, data, momentargs);
    m2 = iv(theta2, data, momentargs);
	m = [m1 m2];
endfunction

# The OLS estimator
function m = qml(theta, data, momentargs)
	k = momentargs{1};
	y = data(:,1);
	x = data(:,2:k+1);
	w = data(:,k+1:columns(data));
	e = y - x*theta;
	m = diag(e)*x;
endfunction

# the IV estimator
function m = iv(theta, data, momentargs)
	k = momentargs{1};
	y = data(:,1);
	x = data(:,2:k+1);
	w = data(:,k+2:columns(data));
	e = y - x*theta;
	m = diag(e)*w;
endfunction

# generate data
[y, z, w] = dgp_linear(n, good_inst, rho);
z = [ones(n,1) , z];
[w, junk, junk2] = st_norm(w);
w = [ones(n,1) , w];


# define things for GMM estimation
data = [y, z, w];
moments = "joint";
k = columns(z);
momentargs = {k}; # tells the moment function where regressors stop and instruments start
b_ols = ols(y,z); # start values
theta = [b_ols ; b_ols]; # use 2 copies of OLS results for GMM start values
control = {100,2}; # bfgsmin controls

# first round inefficient GMM to get consistent estimate
weight = eye(k+columns(w));
[theta, obj_value, iterations, convergence] = gmm_estimate(theta, data, weight, moments, momentargs, control);

# get the efficient weight matrix for the joint moments using consistent 1st round estimate
m = feval(moments, theta, data, momentargs); # evaluate moment conditions
omega = cov(m);
weight = inv(omega);

# get the 2nd round efficient joint estimates
[theta2, obj_value, iterations, convergence] = gmm_estimate(theta, data, weight, moments, momentargs, control);

vc = gmm_variance(theta2, data, weight, moments, momentargs);
k = rows(vc)/2;
mask = [ones(k,k) zeros(k,k); zeros(k,k) ones(k,k)];
vc2 = vc.*mask;

# Standar Hausman test using the masked covariance which ignores covariance of estimators
R = [eye(k) ,  -eye(k)];
H = theta2' * R' * inverse(R*vc2*R') * R * theta2;
q = 1; # degrees of freedom (only 1 since model is linear and the OLS and IV estimators share a moment condition)
pvalue = 1 - chi2cdf(H,q);
printf("\nStandard Hausman test value: %f    P-value: %f\n", H, pvalue);

# Modified Hausman test using the correct covariance and overall efficient weight matrix (H2 test)
R = [eye(k) ,  -eye(k)];
H = theta2' * R' * inverse(R*vc*R') * R * theta2;
q = 1; # degrees of freedom (only 1 since model is linear and the OLS and IV estimators share a moment condition)
pvalue = 1 - chi2cdf(H,q);
printf("\nModified Hausman test value: %f    P-value: %f\n", H, pvalue);
