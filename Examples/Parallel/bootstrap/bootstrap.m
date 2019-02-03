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
 

# bootstrap.m: Obtains reps bootstrap evaluations of the function f
#
# usage: output = bootstrap(f, f_args, reps, nslaves)
# 
# inputs:
#  f: name of function to bootstrap. Should have f_args (cell)
#			as ONLY argument
#  f_args: arguments to the function. (1,1) element is the data,
#     an nxk matrix, where the n rows hold IID observation vectors
#  reps: the number of replications desired
#  nslaves: the number of slave nodes, 0 for serial evaluation on 1 computer
#
# f should return a row vector of output from feval(f, f_args)

function output = bootstrap(f, f_args, reps, nslaves)
	
	# determine size of output vector from f, and create container for results
	output = feval(f, f_args);
	n_returns = size(output,2);
	output = zeros(reps, n_returns);

	# check if doing this parallel or serial
	global PARALLEL = 0; # default is serial
  
  	if (nargin == 4)
		if nslaves > 0
			global NSLAVES NEWORLD NSLAVES TAG;
		  	PARALLEL = 1;
	    		NSLAVES = nslaves;
		endif
	endif

	if !PARALLEL # ordinary serial version
		data = f_args{1};
		for i = 1:reps
			new_data = bootstrap_resample_iid(data); # resample
			f_args{1} = new_data; # evaluate f with resample data
			output(i,:) = feval(f, f_args);
		endfor
	else # parallel version
		LAM_Init(nslaves);

		# The command that the slave nodes will execute
  		cmd=[	'contrib = bootstrap_nodes(f, f_args, n_returns, nn); ',...	
      	 		'MPI_Send(contrib,0,TAG,NEWORLD);'];	

    		nn = floor(reps/(NSLAVES + 1));	# How many reps per slave? Rest is for master

		NumCmds_Send({'f', 'f_args', 'n_returns', 'nn', 'cmd'}, {f, f_args, n_returns, nn, cmd}); # Send all data to all nodes

		# run command locally for last block (slaves are still busy)
		n_master = reps - NSLAVES*nn; # how many to do?
		contrib = bootstrap_nodes(f, f_args, n_returns, n_master);
  		output(reps - n_master + 1:reps,:) = contrib;

 		# collect slaves' results
		contrib = zeros(nn, n_returns);
  		for i = 1:NSLAVES
	  		MPI_Recv(contrib,i,TAG,NEWORLD);
			startblock = i*nn - nn + 1;
    			endblock = i*nn;
  			output(startblock:endblock,:) = contrib;  	
  		endfor

		# clean up after parallel
		if PARALLEL LAM_Finalize; endif

	endif
	
endfunction
