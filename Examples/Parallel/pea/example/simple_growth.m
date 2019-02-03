# Copyright (C) 2005, 2007  Michael Creel <michael.creel@uab.es>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# simple_growth.m: example model for solution by pea.m. This is illustrative
# of the i/o format the model needs to have. This simulates the simple stochastic
# growth model. There is a C++ version (simple_growth.cc) that does the same thing,
# and which is approximately 7-10 X faster
function [expectations_data, hit] = simple_growth(exp_params, model_params, randseed, pea_args, up_bound, low_bound)
	# model parameters
	alpha   = model_params(1,:);  # Capital share
	delta   = model_params(2,:);   # Discount factor
	gam     = model_params(3,:);   # Risk aversion parameter
	d       = model_params(4,:);   # Depreciation rate
	sigma   = model_params(5,:);   # Standard deviation for log noise
	rho     = model_params(6,:);  # Persistence of log technology shock

	# pea controls
	T = pea_args{1};
	burnin = pea_args{2};
	nodes = pea_args{6} + 1; # total nodes

	T = ceil(T / nodes)+1; # number of periods to do on each node (round up if uneven)
	randn("seed",randseed); # set the seed to the fixed value to avoid chatter

	# make containers and initialize
	epsi = randn(T + burnin,1)*sigma;
	shocks = zeros(T + burnin,1);
	k = zeros(burnin+T+1,1);
 	ks = ((1-delta+delta*d) / (alpha*delta) )^(1/(alpha-1));
	k(1,:) = ks;
	c = zeros(burnin+T,1);
	hit = 0;

	# simulate, with burnin
	for t = 1:(burnin+T)
		if (t==1)
			shocks(t) = 1;
		else
			shocks(t,:) = (shocks(t-1) ^ rho) * exp(epsi(t));
		endif
		uprime = exp(exp_params(1) + exp_params(2)*log(k(t)) + exp_params(3)*log(shocks(t)));
		c(t) = (delta*uprime )^(-1/gam);
		k(t+1) = k(t)^alpha*shocks(t) - c(t) + (1-d)*k(t);
		if k(t+1) > up_bound
    		k(t+1) = k(1); hit = 1;
		elseif k(t+1) < low_bound
    		k(t+1) = k(1); hit = 1;
		endif;
		c(t) = k(t)^alpha*shocks(t)  + (1-d)*k(t)-k(t+1);
	endfor;

	# chop off the burnin period
	k = k(burnin + 1:rows(k));
	c = c(burnin + 1: rows(c));
	shocks = shocks(burnin + 1:rows(shocks));
	# data for model of expectations
	e = zeros(T-1,1);
	for t = 1:T-1
		e(t) = c(t+1)^(-gam)*(1-d+k(t+1)^(alpha-1)*alpha*shocks(t+1));
	endfor;
	k = k(1:T-1);
	shocks = shocks(1:T-1);
	x = [ones(T-1,1) log(k) log(shocks)];
	expectations_data = [e x];
endfunction
