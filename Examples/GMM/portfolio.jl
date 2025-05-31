#
# This is an example of estimation of a portfolio model
# as in Hansen and Singleton, 1982. The data comes from
# Tauchen, JBES, 1986 (available by ftp).
# The notation follows my lecture notes. In particular
# gamma is the coefficient of RRA. Try out different sets
# of intruments and different lag lengths. You will see
# that beta is stable across  models, but that gamma is not.
#
# Michael Creel, Oct. 24, 2000
# revised 7 Nov. 2002
# translated to Octave 9/9/2003
# translated to Julia 07 Aug. 2017
# michael.creel@uab.es


# The main thing you have to to is define the moments for estimation,
# as in the function that follows immediately, and ensure that the
# efficient weight matrix is estimated appropriately. The rest
# is already programmed.


# the following function defines the moment conditions, given the parameter
# vector and the instruments. It returns the individual contributions
# to allow estimation of the covariance matrix of the moments.
using DelimitedFiles, Statistics, Econometrics
function portfolio_moments(theta, data)
	# parameters
	beta = theta[1]
	gam = theta[2]
	# data items
	c = data[:,1]
	r = data[:,2]
    inst = data[:,3:end]
	#  form error function
	# note that c = c_t / c_t-1 (for stationarity) was done in data preparation
	e = 1.0 .- beta*(1.0 .+ r) .* (c .^ (-gam))
	# cross with instruments
	m = e.*inst
end

## read the data and form variables
cd(@__DIR__)
data = readdlm("../Data/tauchen.data")
c = data[:,1]
p = data[:,2]
d = data[:,3]
# form net return and stationary consumption
r = (p + d) ./ lag(p,1) .- 1.0
c = c ./ lag(c,1); # ensure stationarity

## choose the instruments: note - the code is off by one w.r.t. the notes
inst = [lag(c,1) lag(r,1)]
maxlag = 5
inst = [inst lags(inst, maxlag)]
inst = inst[maxlag+2:end,:]  # drop first two obs, due to two lags of c
inst, junk, junk = stnorm(inst) # for numeric stability
inst = [ones(size(inst,1),1) inst]  # add a constant
# collect the data for the moments
data = [c[maxlag+2:end] r[maxlag+2:end] inst]

## do estimation
moments = θ -> portfolio_moments(θ, data)
weight = 1.0
names = ["beta","gamma"]
# initial consistent estimate
θstart  = [0.95, 2.]
θhat, obj_value, D, ms, convergence = gmm(moments, θstart , weight)
# second step with efficient weight matrix
Ω  = cov(ms)
weight = inv(Ω)
title = "Two step GMM estimation of portfolio model"
results = gmmresults(moments, θhat, weight, title, names)

## CUE GMM
title = "CUE GMM estimation of portfolio model"
gmmresults(moments, θstart , "", title, names)
nothing
