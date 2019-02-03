## Copyright (C) 2013 Michael Creel
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
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

## EstimateBVAR

## Author: Michael Creel <michael@yosemite>
## Created: 2013-12-18

function [ ret ] = EstimateBVAR ()
	close all;
	load rbcdata.m; % columns are cons, inv, hours
	data = rbcdata;

	% demean and normalize the data, makes homosced VAR more plausible 
	% also, if the variables have marginal variance = 1, we can
	% set the prior precision in relation to this
	data = st_norm(data);
	PP = 0.5; % prior precision: small imposes weakly

	maxlag = 2; % how many lags to use?

	% set up the model
	data = [data lags(data,maxlag)]; % add lags
	data = data(maxlag+1:end,:); % drop rows with missing
	y = data(:,1:3);
	plot(y);
	legend('c','i','n');
	title('the data');
	figure;
	x = data(:,4:end);

	% First plain OLS
	b = ols(y,x);
	printf('plain OLS results\n');
	printf('A1\n');
	disp(b(1:3,:));
	printf('A2\n');
	disp(b(4:6,:));
	fit = x*b;
	e = y-fit;
	ess = sum(e.*e);
	tss = sum(y.*y);
	rsq = 1-ess./tss;
	printf('r-squares OLS fit: %f\n');
	disp(rsq);
	n = rows(y);
	printf('\n\n');
	% now Minnesota priors
	y = [y; PP*eye(3)]; 	% first the restriction on first lag coefficients
	if maxlag > 1 		% now the higher order lags restricted to zero
		for j=2:maxlag
			y = [y; PP*zeros(3,3)];
		endfor
	endif	
	x = [x; PP*eye(3*maxlag)];
	b = ols(y,x);
	printf('#################################\n');
	printf('Minnesota prior results\n');

	printf('A1\n');
	disp(b(1:3,:));
	printf('A2\n');
	disp(b(4:6,:));
	x = x(1:n,:);
	y = y(1:n,:);
	fit = x*b;
	e = y-fit;

	ess = sum(e.*e);
	tss = sum(y.*y);
	rsq = 1-ess./tss;
	printf('r-squares Minnesota fit: %f\n');
	disp(rsq);
	plot(fit);
	legend('c','i','n');
	title('the Minnesota fit');
	figure;

	corr(e)
	
	e = e + [-1 0 1];
	plot(e);
	legend('c','i','n');
	title('residuals (with separation added in), Minnesota fit');

endfunction
