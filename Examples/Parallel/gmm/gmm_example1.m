load data;  # read data

k = 5; # number of regressors
theta = zeros(k,1); # start values
weight = eye(columns(data) - 1 - k); # weight matrix
moments = "Poisson_IV_moments"; # name of function that calculates moments
momentargs = {k}; # additional information needed to calculate moments
control = {100,0};
# gmm estimation: serial
theta_s = gmm_estimate(theta, data, weight, moments, momentargs, control);

# gmm estimation: parallel
nslaves = 1;
theta_p = gmm_estimate(theta, data, weight, moments, momentargs, control, nslaves);

printf("GMM parameter estimates\nMoment conditions: %s\n%d observations\n", moments, rows(data));

# print results
labels = char("serial", "parallel");
prettyprint_c([theta_s theta_p], labels);

