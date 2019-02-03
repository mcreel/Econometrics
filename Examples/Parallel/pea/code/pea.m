# Copyright (C) 2006,2007  Michael Creel <michael.creel@uab.es>
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

# pea: solves a nonlinear rational expectations model using the method of parameterized expectations
# usage:  output = pea(model, model_params, exp_model, exp_params, pea_args)
# inputs:
#   model: (string) name of function that contains the model to be solved
#   model_params: vector that holds the model's parameters
#   exp_model: (string) name of function that models expectations
#   exp_params: vector of parameters for the model of expectations
#   pea_args: 7 element cell array of controls for pea. These are:
#       T: the length of simulation to do
#       burnin: number of periods to pre-pend. These are discarded on all nodes.
#       start_value: initial state variable(s)
#       homotopy: a scalar in (0,1]. Low values give slow update. 1 gives no memory.
#       maxiters: maximum number of iterations (number of re-calculations of model of expectations)
#       nodes: number of compute nodes to use. Set to zero for serial.
#       verbosity: 0,1,2, or 3. Increasing levels of screen output.
# outputs:
#   output: a row vector. The first element is a 0/1 indicator of convergence (1=yes)
#           The other elements are the final values of the parameters of the model of expectations.
#
# Implements a moving bounds PEA algorithm, as described in
# Lilia Maliar and Seguei Maliar "Parameterized Expectations Algorithm
# and the Moving Bounds",  Journal of Business and Economic Statistics,
# 21, January 2003, pp. 88-92. MATLAB code provided by Maliar and Maliar
# was used as a model for this code, but the present code is essentially
# a new implementation of the algorithm.
#
# Details: users must supply two functions: a model to be solved, and a model of
# expectations. The provided example functions simple_growth(.cc or .m) and simple_expectations.m
# illustrate the format that must be used in order for pea.m to solve the model.
function output = pea(model, model_params, exp_model, exp_params, pea_args)
	# pea controls
	T = pea_args{1};
	start_value = pea_args{3};
	homotopy = pea_args{4};
	maxiters = pea_args{5};
	compute_nodes = pea_args{6};
	verbose = pea_args{7};
	control = {4, verbose}; # bfgsmin controls for NLS fitting expectations to data
	global NEWORLD PARALLEL TAG NSLAVES
	PARALLEL = 0;
	# Initialize MPI if needed
	if compute_nodes > 0
		NSLAVES = compute_nodes;
		PARALLEL = 1;
		LAM_Init(compute_nodes);
	endif
	# how many periods for each node?
	nodes = compute_nodes + 1;
	T = ceil(T / nodes); # number of periods to do on each node (round up if uneven)
	# send things out to nodes, if parallel
	if PARALLEL
		# send static items to nodes
		pea_args{1} = T; # split it up if parallel
		cmd = ["randseed=randn('seed');"];
		NumCmds_Send({"model", "model_params", "exp_model", "exp_params", "pea_args", "cmd"},{model, model_params, exp_model, exp_params, pea_args, cmd});
	endif
	# algorithm parameters
	crate = 0.007; # Speed of moving the bound
	criter = 1e-5; # Convergence criterion
	# counters and controls
	iteration  = 1;
	dif = 2e-5; # tolerance for parameter convergence
	hit = 0;
	# MAIN LOOP
	tic;
	converged = 0;
	randseed = randn("seed"); # need to keep it fixed to avoid chatter across iterations
	while (!converged) & (iteration < maxiters)
		# cmd has to be re-set each iteration since NLS uses a different cmd
		cmd = ["[exp_data, hit_node] = feval(model, exp_params, model_params, randseed, pea_args, up_bound, low_bound);",...
			"MPI_Send(hit_node, 0, TAG, NEWORLD);"];
		up_bound = start_value*(2-exp(-crate*iteration)); # Upper bound
		low_bound = start_value*exp(-crate*iteration);     # Lower bound
		hit = 0; # Indicator if bound is hit
		# send parameters of expectations model and new bounds to slaves
		if PARALLEL
			NumCmds_Send({"exp_params", "up_bound", "low_bound", "cmd"}, {exp_params, up_bound, low_bound, cmd});
		endif
		# run first block on master while slaves are busy
		[exp_data, hit] = feval(model, exp_params, model_params, randseed, pea_args, up_bound, low_bound);
		# collect slaves' results
		if PARALLEL
			hit_node = 0.0;
			for i = 1:compute_nodes
				MPI_Recv(hit_node, i, TAG, NEWORLD);
				hit = hit + hit_node;
			endfor
		endif
		# NLS fit to update expectation function's parameters
		# Note that exp_data is local to compute nodes, already in their memory. No need to send it
		# On the other hand, exp_params needs to be sync'ed across nodes. This is done by pea_obj.m
		new_params = bfgsmin("pea_obj", {exp_params, exp_model, exp_data, compute_nodes}, control);
		dif = norm(exp_params-new_params); # criterion value
		converged = ((dif <= criter) & (hit==0));  # convergence required negligable param change and no bounds hit
		exp_params = homotopy*new_params + (1-homotopy)*exp_params; # update the coefficients
		iteration = iteration + 1; # iteration couter
		hit = hit > 0;
	endwhile;
	t = toc;
	printf("Time %f seconds for %d iterations on %d nodes\n", t, iteration, compute_nodes+1);
	if (compute_nodes > 0) LAM_Finalize; endif
	output = [converged exp_params'];
endfunction
