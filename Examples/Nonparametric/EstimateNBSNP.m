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
x_orig = x;
# choose one of the next 2 lines to see the effects of scaling the data
[x scale] = scale_data(x); # scale the data
data = [y x ];
names = char("constant","pub. ins.","priv. ins.", "sex", "age","edu","inc","gam1", "gam2", "lnalpha");
model = "NegBinSNP";
nb_type = 1;
modelargs = {nb_type};
theta = zeros(columns(x) + 3,1);
mletitle = "NegBin SNP model, MEPS full data set";
# SA controls
  ub = 2*ones(rows(theta),1);
  lb = -ub;
  nt = 5;
  ns = 2;
  rt = 0.75; # careful - this is too low for many problems
  maxevals = 1e4;
  neps = 5;
  functol = 1e-10;
  paramtol = 1e-3;
  verbosity = 2;
  minarg = 1;
  control = { lb, ub, nt, ns, rt, maxevals, neps, functol, paramtol, verbosity, 1};
theta = mle_estimate(theta, data, model, modelargs, control);


# BFGS controls
control = {-1,2};
theta = mle_results(theta, data, model, modelargs, names, mletitle, scale, control);



