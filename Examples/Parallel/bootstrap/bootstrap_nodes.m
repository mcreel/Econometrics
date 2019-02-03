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


# this does a portion of the replications on each slave
function contrib = bootstrap_nodes(f, f_args, n_returns, nn)
	contrib = zeros(nn, n_returns); # container for results
	data = f_args{1}; # get the real data for resample
	for i = 1:nn
		new_data = bootstrap_resample_iid(data); # resample		
		f_args{1} = new_data; # evaluate f with resample data
		contrib(i,:) = feval(f, f_args);
	endfor
endfunction		
