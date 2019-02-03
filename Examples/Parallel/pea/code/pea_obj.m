# Copyright (C) 2005,2007  Michael Creel <michael.creel@uab.es>
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

# pea_obj: for internal usage by pea.m

function [obj_value, score] = pea_obj(exp_params, exp_model, exp_data, compute_nodes)

	if compute_nodes > 0
		global NEWORLD NSLAVES TAG PARALLEL

		# The command that the compute nodes will execute
    		cmd=[	'contrib = feval(exp_model, exp_params, exp_data);',...
			'MPI_Send(contrib, 0, TAG, NEWORLD);'
			];

		# send items to slaves
		NumCmds_Send({"exp_params", "cmd"}, {exp_params, cmd});

		# evaluate last block on master while compute nodes are busy
	  	running_sum = feval(exp_model, exp_params, exp_data);
		# collect nodes' results
		contrib = zeros(size(running_sum));
	  	for i = 1:NSLAVES
			MPI_Recv(contrib, i , TAG, NEWORLD);
			running_sum = running_sum + contrib;
		endfor
	else # serial version
 		running_sum = feval(exp_model, exp_params, exp_data);
  	endif

	# define pieces
	nodes = compute_nodes + 1; # number of nodes
	T = rows(exp_data) * nodes; # total number of data points
	running_sum = running_sum / T; # put into averages, for better numeric stability
	obj_value = running_sum(1,1);
	score = running_sum(2:rows(running_sum),:);

	# let's bullet-proof this in case the model goes nuts
	if (((abs(obj_value) == Inf)) || (isnan(obj_value)))
		obj_value = realmax/2;
		score = "na";
	endif
endfunction
