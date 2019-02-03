## Copyright (C) 2010 Michael Creel
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## AR1Errors

## Author: Michael Creel <michael@yosemite>
## Created: 2010-03-11


# this script does a Monte Carlo that shows efficiency of GLS compared to OLS
# when there are AR1 errors

1;
# the function to be Monte Carlo'ed
function b = wrapper(args)
	n = args{1};    # sample size
	rho = args{2};  # AR1 parameter
	x = randn(n,1); # an exogenous regressor
	# make the AR1 errors
	e = zeros(n,1);
	e(1,:) = randn(1,1)/sqrt(1-rho^2);
	for t = 2:n
	  e(t,1) = rho*e(t-1,:) + randn(1,1); # the white noise
	endfor
	# the dep var
	y = 1 + x + e;
	# OLS
	x = [ones(n,1) x];
	b_ols = ols(y,x)';
	# estimate rho
	e = y - x*b_ols';
	rhohat = ols(e(2:n,:),e(1:n-1,:));
	# FGLS 
	ystar = y - rhohat*lag(y,1);
	xstar = x - rhohat*lag(x,1);
	ystar(1,:) = sqrt(1-rhohat^2)*y(1,:);
	xstar(1,:) = sqrt(1-rhohat^2)*x(1,:);
	b_gls = ols(ystar, xstar)';
	b_ols = b_ols - [1 1]; # subtract true values, so mean should be approx. zero if consistent
	b_gls = b_gls - [1 1]; # subtract true values, so mean should be approx. zero if consistent
	b = [b_ols b_gls];
endfunction

# do the Monte Carlo
n = 30;
rho = 0.9;
args = {n, rho};
outfile = 'AR1Errors.out';
n_pooled = 100;
verbose = true;
reps = 1000;
system('rm AR1Errors.out'); # start from scratch each time
montecarlo('wrapper', args, reps, outfile, n_pooled, false, verbose);


# analyze results
load AR1Errors.out;
results = AR1Errors(:,3:6); # drop the timing and node info added by montecarlo.m
clear AR1Errors;
close all

hist(results(:,2),30);
#print('AR1ErrorsOLS.png', '-dpng');
figure;
hist(results(:,4),30);
#print('AR1ErrorsGLS.png', '-dpng');
dstats(results);
		

