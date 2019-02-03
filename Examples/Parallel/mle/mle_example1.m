load poisson_data;
model = "Poisson"; 	# name of function that calculates loglikelihood

theta = zeros(columns(poisson_data)-1,1); # starting values
control = {100, 0};

# Estimation: serial
t = cputime();
theta_s = mle_estimate(theta, poisson_data, model, "", control);
cputime() - t
# Estimation parallel
nslaves = 1;
t = cputime();
theta_p = mle_estimate(theta, poisson_data, model, "", control, nslaves);
cputime() - t
printf("MLE parameter estimates: %s model, %d observations\n", model, rows(poisson_data));
labels = char("serial", "parallel");
prettyprint_c([theta_s theta_p], labels);
