# Copyright (C) 2007  Michael Creel <michael.creel@uab.es>
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

# parallel_performance: example code that looks at parallel performance of pea.m
# usage: examine it and run it! (you need to have MPITB installed and a LAM running)

# model
model = "simple_growth";
alpha   = 0.33;  # Capital share
delta   = .95;   # Discount factor
gam     = 1.0;   # Risk aversion parameter
d       = .02;   # Depreciation rate
sigma   = .01;   # Standard deviation for log noise
rho     = 0.95;  # Persistence of log technology shock
model_params = [alpha; delta; gam; d; sigma; rho];

# model of expectations
exp_model = "simple_expectations";
exp_params = zeros(3,1);;

# pea_args
T = 200000; # length of simulation
burnin = 100; # initial periods to drop

# initial value of of bounds, here I'm using steady state for k
start_value = ((1-delta+delta*d) / (alpha*delta) )^(1/(alpha-1));
homotopy = 1.0; # set this to a number in (0,1]. 1 is no memory, closer to 0 is slow update.
maxiters = 30; # limit to iterations (useful for timings)
nodes = 0; # number of compute nodes to use (0 for serial)
verbosity = 0; # observe intermediate output (0,1,2 or 3. 3 is maximum detail)
pea_args = {T, burnin, start_value, homotopy, maxiters, nodes, verbosity};

# loop over several cluster sizes
printf("Sample size: %d burnin: %d  maxiters %d\n", T, burnin, maxiters);
for nodes = 0:1
	pea_args{6} = nodes;
	pea(model, model_params, exp_model, exp_params, pea_args);
endfor
