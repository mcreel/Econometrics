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

## load the data
using Econometrics, DelimitedFiles
which_dep = 1
if (which_dep == 1) dep = "OBDV"   end
if (which_dep == 2) dep = "OPV"    end
if (which_dep == 3) dep = "IPV"    end
if (which_dep == 4) dep = "ERV"    end
if (which_dep == 5) dep = "DV"     end
if (which_dep == 6) dep = "PRESCR" end
cd(@__DIR__)
data = readdlm("../Data/meps1996.data")

## create the regressors and dep var
y = data[:,which_dep]
x = data[:,7:12]
n = size(x,1)
x = [ones(n,1) x]

## to check effects of poor scaling, try commenting the next line
# Note: mle.jl now uses automatic differentiation by default, which
# is not sensitive to scaling, in this case. So, you will not see a difference
# for the OBDV data. Perhaps, for others, there will be a difference, haven't checked
x[:,end] = x[:,end]/1000.0 # put income in thousands

## do ML
names = ["constant", "pub. ins.","priv. ins.", "sex", "age","edu","inc"]
title = "Poisson model, "*dep* ",  MEPS 1996 full data set"
model = theta -> poisson(theta, y, x)
θstart = zeros(size(x,2)) # start values for estimation
# NOTE: vc=1 uses default sandwich covariance estimator.
# try adding the option vc=2 for "Hessian" or vc=3 for OPG
θhat, objvalue, V, converged = mleresults(model, θstart, title, names, vc=1);
