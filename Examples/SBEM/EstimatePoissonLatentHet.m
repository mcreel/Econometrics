compute_nodes = 1;

# This does simple Poisson estimation using the meps1996.data

# The MEPS data

# Define dep and expl vbls
# 
# 	The dep. vbls, in corresponding column
# 	Office based doctor visits	1
# 	Outpatient doctor visits	2
# 	Emergency room visits		3
# 	Inpatient visits		4
# 	Dental visits			5
# 	Prescriptions			6
# 	

1;


which_dep = 1;
if (which_dep == 1) printf("\nOBDV\n"); endif
if (which_dep == 2) printf("\nOPV");endif
if (which_dep == 3) printf("\nIPV");endif
if (which_dep == 4) printf("\nERV");endif
if (which_dep == 5) printf("\nDV");endif
if (which_dep == 6) printf("\nPRESCR");endif

load meps1996.data;
y = data(:,which_dep);
x = data(:,7:12);
n = rows(x);
x = [ones(n,1) x];
[x scale] = scale_data(x); # scale the data
k = columns(x);
S = 1000; # number of simulations

# start values (these are good values obtained from a previous run, to speed things up)
theta = [
0.331422;
0.358579;
0.304063;
0.207217;
0.136424;
0.142148;
-0.032036;
0.146271];
model = "PoissonLatentHet";
otherargs = {k};
data = [y x randn(n, S)];
names = char("constant","pub. ins.","priv. ins.", "sex", "age","edu","inc", "lnalpha");
title = "Poisson Latent Heterogeneity model, SML estimation, MEPS 1996 full data set";

# samin controls - these are very extreme, since very good start values are used
# in normal practice, one would need to use wide bounds and slower cooling
ub = theta + 0.1;
lb = theta - 0.1;
nt = 2;
ns = 1;
rt = .25;
maxevals = 1e6;
neps = 3;
functol = 1e-7;
paramtol = 1e-4;
verbosity = 2;
minarg = 1;
control = { lb, ub, nt, ns, rt, maxevals, neps, functol, paramtol, verbosity, 1};
theta = mle_estimate(theta, data, model, otherargs, control);

# iters set to zero to avoid attempting BFGS min of non-differentiable fn.
control = {0};
mle_results(theta, data, model, otherargs, names, title, scale, control);

