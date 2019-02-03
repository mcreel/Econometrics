function results = wrapper(args)
	model = args{1};
	theta = args{2};
	n = args{3};
	simlength = args{4};
	burnin = args{5};
	vc_reps = args{6};
	from_param_space = args{7}; # use true value, or a draw from parameter space?

	# move thetastart up here so it can be used as true value
	ub = [2; 0.999; 1; 1];
	lb = [0; 0; 0; -2];

	thetastart = rand(size(theta)).*(ub-lb) + lb;
	
	if from_param_space theta = thetastart; endif

	# make Monte Carlo data
	modelargs = {"", n, burnin, true};
	[data, L, K, snm_draws] = feval(model, theta, modelargs);

	# define a few things
	M = columns(data) - L - K;

	ub = [ub; 3*ones(K,1)];
	lb = [lb; 0.1*ones(K,1)];
	thetastart = rand(size(ub)).*(ub-lb) + lb;
	
	# prepare data, with trimming using kernel density of condvars
	condvars = data(:,L+1:L+K);
	d = kernel_density(condvars, condvars);
	keep = d > quantile(d, 0.02);
	data = data(keep,:);
	endogs = data(:,1:L);
	condvars = data(:,L+1:L+K);
	instruments = data(:,L+K+1:L+K+M);
	n = rows(data);
	
	# prewhiten and divide by ww
	[endogs, condvars, P] = snm_dataprep(endogs, condvars);

	data = [endogs condvars instruments];

	# define momentargs for GMM, including random draws for long simulation
	modelargs = {"", simlength, burnin, false};
	[junk, junk, junk, snm_draws] = feval(model, theta, modelargs);
	modelargs{1} = snm_draws;
	momentargs = {model, modelargs, L, K, M, P, true, false};

	# samin controls
	nt = 3;
	ns = 3;
	rt = .25;
	maxevals = 1e4;
	neps = 3;
	functol = 1e-6;
	paramtol = 1e-3;
	verbosity = 2;
	minarg = 1;
	sacontrol = { lb, ub, nt, ns, rt, maxevals, neps, functol, paramtol, verbosity, 1};	

	# first round consistent estimator
	weight = 1; # weight matrix


#f = @(x) gmm_obj(x', data, weight, "snm_moments", momentargs, 0); 
# 
#ctl.XVmin = lb';
#ctl.XVmax = ub';
#ctl.tol = 1e-1;
#ctl.refresh = 1;
#ctl.constr = 1;
#ctl.maxnfe = 500;
#[thetahat, obj_value, junk, conv] = de_min (f, ctl);

	[thetahat, obj_value, conv1] = gmm_estimate(thetastart, data, weight, "snm_moments", momentargs, sacontrol, 0);
	thetahat1 = parameterize(thetahat)

	# get feasible efficient weight matrix using first round
# 	verbose = false; # output from cov matrix routine?
# 	omega = snm_covariances(thetahat, data, momentargs, vc_reps, verbose);
#	weight = inv(diag(diag(omega)));
# 	[thetahat, obj_value, conv2] = gmm_estimate(thetahat, data, weight, "snm_moments", momentargs, sacontrol, 0);
# 	thetahat2 = parameterize(thetahat)

	 
	results = [theta' thetahat1(1:4,:)' thetahat1(1:4,:)' obj_value];
endfunction
