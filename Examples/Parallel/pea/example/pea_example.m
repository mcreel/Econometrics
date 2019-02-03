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

# pea_example: example code that shows how to use pea.m
# usage: examine it and run it!

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
maxiters = 1e10; # limit to iterations (useful for timings)
nodes = 0; # number of compute nodes to use (0 for serial)
verbosity = 2; # observe intermediate output (0,1,2 or 3. 3 is maximum detail)
pea_args = {T, burnin, start_value, homotopy, maxiters, nodes, verbosity};

# run it
output = pea(model, model_params, exp_model, exp_params, pea_args);

# make some plots
exp_params = output(:,2:columns(output))';
pea_args{1} = 2003;
data = feval(model, exp_params, model_params, randn("seed"), pea_args, Inf, 0);
k = exp(data(:,3));
shocks = exp(data(:,4));
T = rows(k);
c = zeros(T-1,1);
for t = 1:T-1
	c(t) = k(t)^alpha*shocks(t,1) + (1-d)*k(t)-k(t+1);
endfor

# trim to get rid of missings
T = T - 1;
investment = k - (1-d)*lag(k,1);
investment = investment(2:T,:);
consumption = c(2:T,:);
k = k(2:T,:);
shocks = shocks(2:T,:);
output = shocks.*k.^alpha;
T = T - 1;

time=(1:1:T);
subplot(2,2,1);
title('Time series solution');
xlabel('t');

plot(time,k);
ylabel('Capital');

subplot(2,2,2);
plot(time,consumption);
ylabel('Consumption');

subplot(2,2,3);
plot(time,investment);
ylabel('Investment');

subplot(2,2,4);
plot(time, output);
ylabel('Output');

