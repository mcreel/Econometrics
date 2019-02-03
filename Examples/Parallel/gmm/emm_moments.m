# Copyright (C) 2006  Michael Creel <michael.creel@uab.es>
# under the terms of the GNU General Public License.
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


# defines EMM moments.
function scores = emm_moments(theta, data, momentargs)

	k = momentargs{1};
	dgp = momentargs{2}; # the data generating process (DGP)
	dgpargs = momentargs{3}; # its arguments (cell array)
	sg = momentargs{4}; # the score generator
	sgargs = momentargs{5}; # SG arguments (cell array)
	phi = momentargs{6}; # QML estimate of SG parameter

	y = data(:,1);
	x = data(:,2:k+1);
	
	# random draws passed in data to ensure fixed over estimation
	rand_draws = data(:,k+2:columns(data));

	n = rows(y);
	
	scores = zeros(n,rows(phi));
	reps = columns(rand_draws);
	
	for i = 1:reps
		e = rand_draws(:,i);
		y = feval(dgp, theta, x, e, dgpargs); # simulated data
		sgdata = [y x]; # pass simulated data to SG
		scores = scores + numgradient(sg, {phi, sgdata, sgargs}); # gradient of SG
	endfor

	scores = scores / reps; # average over number of simulations

endfunction
 
