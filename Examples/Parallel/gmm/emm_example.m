# This shows how EMM estimation may be done.
# The DGP is a probit model, and we estimate by EMM using a logit
# model as the score generator. The simulated binary responses are
# smoothed using an approximation to a step function. The bigger that
# you set dgpargs, the better the approximation to a step function.


n = 1000;
k = 3;
reps = 10;  # simulations in EMM

dgp = "ProbitDGP"; # the DGP
dgpargs = 4; # smoothing parameter, to make objective fn differentiable
sg = "Logit"; # the score generator
sgargs = {1}; # just a place holder in this example

theta = ones(k,1);

x = [ones(n,1) randn(n,k-1)];
y = ProbitDGP(theta, x);
data = [y x];

# QML estimation of score generator
control = {100,2,1,1}; # BFGS controls
printf("fitting the score generator\n");
%phi = bfgsmin("mle_obj", {theta, data, "Logit",''}, control);
phi = mle_estimate(theta, data, 'Logit','', control);

moments = "emm_moments"; # define moments function for GMM

# arguments for emm_moments
momentargs = {k, dgp, dgpargs, sg, sgargs, phi};

# arguments for emm_estimate
theta = phi; # start values
weight = eye(rows(phi)); # weight matrix
rand_draws = randn(n,reps); # fixed over iterations to avoid "chattering"
data = [data rand_draws];

control = {100,2,1,1}; # BFGS controls
# first round consistent estimator
[theta, obj_value, convergence, iters] = gmm_estimate(theta, data, weight, moments, momentargs, control);
m = emm_moments(theta, data, momentargs);

# now use efficient weights
weight = inverse(cov(m));
names = char("p1","p2","p3");
title = "EMM example";
gmm_results(theta, data, weight, moments, momentargs, names, title, "", control);

