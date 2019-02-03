% this is to solve a model many times using different parameter values
% drawn from a prior distribution. This runs on a cluster. The frontend
% node (MPI rank 0) does not participate in the solution, it only collects
% the results and writes to disk. You should run this using something like
% mpirun -np X --hostfile /foo/bar octave -q --eval make_simdata
% where X is the total number of cores in the cluster, PLUS ONE. That way,
% all cores will be used in computations, and one will have the added duty
% of serving as the frontend

design=true;
verbose=true;

UseDynare;
if not(MPI_Initialized) MPI_Init; end
CW = MPI_Comm_Load("NEWORLD");
node = MPI_Comm_rank(CW);
nodes = MPI_Comm_size(CW);
mytag = 48;


if design
	outfile = "simdata.design";
	reps = 1000;  % total replications
	n_pooled = 10;  % number of runs accumulated before sent from nodes to frontend
else
	outfile = "simdata.paramspace";
	reps = 4e4;  % total replications
	n_pooled = reps/(nodes-1);  % number of runs accumulated before sent from nodes to frontend
endif


parameters; % load information from common, to sync with ASBIL
model_params0 = lb_param_ub(:,2);
lb = lb_param_ub(:,1);
ub = lb_param_ub(:,3);
lb_ub = [lb ub];
if !node
	save lb_ub lb_ub;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% you do not need to alter anything below this line %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
first_time = true;
model_params = model_params0;
% break into pieces
beta  = model_params(1,:);
gam   = model_params(2,:);
rho1   = model_params(3,:);
sigma1 = model_params(4,:);
rho2   = model_params(5,:);
sigma2 = model_params(6,:);
nss   = model_params(7,:);

alpha = 0.33;
delta = 0.025;

% the psi implied by other parameters
c1 = ((1/beta + delta - 1)/alpha)^(1/(1-alpha));
kss = nss/c1;
css = kss * (c1^(1-alpha) - delta);
c2 = (css)^(-gam/alpha);
psi = (1-alpha)*((c1/c2)^(-alpha));

save parameterfile  beta gam rho1 sigma1 rho2 sigma2 nss;

% get the RNG state on this node, to re-establish separate
% states on nodes after running Dynare, which synchronizes them
RNGstate = rand('state');

% solve model once on each node, to get Dynare structures ready
% for simulations
if node==1
	dynare SimpleModel1 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==2
	dynare SimpleModel2 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==3
	dynare SimpleModel3 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==4
	dynare SimpleModel4 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==5
	dynare SimpleModel5 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==6
	dynare SimpleModel6 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==7
	dynare SimpleModel7 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==8
	dynare SimpleModel8 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==9
	dynare SimpleModel9 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==10
	dynare SimpleModel10 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==11
	dynare SimpleModel11 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==12
	dynare SimpleModel12 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==13
	dynare SimpleModel13 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==14
	dynare SimpleModel14 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==15
	dynare SimpleModel15 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==16
	dynare SimpleModel16 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==17
	dynare SimpleModel17 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==18
	dynare SimpleModel18 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==19
	dynare SimpleModel19 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==20
	dynare SimpleModel20 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==21
	dynare SimpleModel21 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==22
	dynare SimpleModel22 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==23
	dynare SimpleModel23 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==24
	dynare SimpleModel24 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==25
	dynare SimpleModel25 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==26
	dynare SimpleModel26 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==27
	dynare SimpleModel27 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==28
	dynare SimpleModel28 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==29
	dynare SimpleModel29 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==30
	dynare SimpleModel30 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==31
	dynare SimpleModel31 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
elseif node==32
	dynare SimpleModel32 noclearall;
	ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
	set_dynare_seed(ss);
end

pernode = reps/(nodes-1);
if node
	for i = 1:pernode/n_pooled
        % break it up to do intermediate writes
		for j = 1:n_pooled
			if design
				model_params = model_params0;
			else
				model_params = rand(size(ub)).*(ub-lb) + lb;
			end
			% break into pieces
            alpha = 0.33;
            delta = 0.025;
			beta  = model_params(1,:);
			gam   = model_params(2,:);
			rho1   = model_params(3,:);
			sigma1 = model_params(4,:);
			rho2   = model_params(5,:);
			sigma2 = model_params(6,:);
			nss   = model_params(7,:);

			% the psi implied by other parameters
			c1 = ((1/beta + delta - 1)/alpha)^(1/(1-alpha));
			kss = nss/c1;
			css = kss * (c1^(1-alpha) - delta);
			c2 = (css)^(-gam/alpha);
			psi = (1-alpha)*((c1/c2)^(-alpha));
			
			% pass values to Dynare
			set_param_value('betta', beta);
			set_param_value('gam', gam);
			set_param_value('rho1', rho1);
			set_param_value('sigma1', sigma1);
			set_param_value('rho2', rho2);
			set_param_value('sigma2', sigma2);
			set_param_value('nss', nss);

			% do the simulation at the param values
			info = stoch_simul(var_list_);
			% get a simulation of length 160 and compute aux. statistic
			data = [y c n r w];
			data = data(101:260,:);
			contrib = vec(data)';
			
			if (i==1) && (j==1) contribs = zeros(n_pooled, columns(contrib)); end
			contribs(j,:) = contrib;
		end
		MPI_Send(contribs, 0, mytag, CW);
	end
else % frontend
	received = 0;
	done = false;
	while received < reps
		% retrieve results from compute nodes
		%pause(0.01);
		for i = 1:nodes-1
			% compute nodes have results yet?
			ready = false;
			ready = MPI_Iprobe(i, mytag, CW); % check if message pending
			if ready
				% get it if it's there
				contribs = MPI_Recv(i, mytag, CW);
				need = reps - received;
				received = received + n_pooled;
				% truncate?
				if n_pooled  >= need
						contribs = contribs(1:need,:);
						done = true;
				end
				FN = fopen (outfile, "a");
				if (FN < 0) error ("make_simdata: couldn't open output file %s", outfile); endif
				for j = 1:rows(contribs)
					fprintf(FN, "%f ", contribs(j,:));
					fprintf(FN, "\n");
				endfor
				fclose(FN);
				system('sync');
				if verbose
					printf("\nContribution received from node%d.  Received so far: %d\n", i, received);
				end
			end
		end
	end
end
MPI_Barrier(CW);
if not(MPI_Finalized) MPI_Finalize; end
