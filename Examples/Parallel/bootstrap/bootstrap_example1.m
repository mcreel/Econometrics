# generates heteroscedastic data and shows that the bootstrap
# st. errors for coefficients are like the Huber-White
# het. consistent standard errors, while the OLS st. errors
# are biased.

# bootstraps the OLS estimator using heteroscedastic data
reps = 2000;
n = 100;
k = 3;

# generate data
x = [ones(n,1) randn(n,k-1)];
e = x(:,k) .* randn(n,1); # heteroscedastic errors
y = x*zeros(k,1) + e;
f = "myols";
f_args = {[y x]}; # we want a cell, remember?

# get bootstrap standard errors
compute_nodes = 0; # this can be done in parallel, but serial by default
output = bootstrap(f, f_args, reps, compute_nodes);

disp("Bootstrap standard errors");
disp(std(output))

# compare to OLS and OLS with consistent varcov
names = char("param1", "param2");
# OLS with standard var-cov estimator
mc_ols(y,x, names, 0, 1);
# OLS with Huber-White var-cov estimator
mc_ols(y,x, names, 0, 0);


