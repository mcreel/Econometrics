# this calculates asymptotic and bootstrap p values for the test
# H0: slope = zero
# The dgp is probit, which is to say, y = 1(0+slope*x > e) where e~N(0,1)

# inputs (in args, a cell):
# n: sample size
# bs_reps: number of bootstrap simulations
# slope: true slope. When this is 0, we're checking size, otherwise, power against H0: slope = 0
# iter_limit: number of bfgs iters allowed for convergence (failures are noted)

# outputs:
# restarts: number of times convergence not obtained within iters_limit
# slope: true slope coefficient of probit dgp
# asymp. p-value
# D-M bootstrap p-value

function results = bs_test(args)

	n = args{1}; # sample size
	bs_reps = args{2}; # number of bootstraps to use
	slope = args{3}; # we're doing size when this is zero, power against H0: slope = 0 otherwize
	iter_limit = args{4};

	# we'll put a limit on iters, and only use trials that converge
	# within the limit
	model = "Probit";
	modelargs = "";

	conv = 0;
	restarts = 0;
	control = {iter_limit, 0};

	while conv != 1
		x = [ones(n,1) -1 + 2*rand(n,1)];
		theta = [0; slope];
		y = ProbitDGP(theta, x);
		data = [y x];
		[thetahat, junk, conv]  = mle_estimate(theta, data, model, modelargs, control);
		if conv != 1 restarts = restarts + 1; endif # tally up convergence failures
	endwhile
	varcov = mle_variance(thetahat, data, model, "");
	se = sqrt(varcov(2,2));
	# power against h0: slope = 0 (this is size when the slope is zero, power otherwise)
	p1 = 2 - 2*normal_cdf(abs(thetahat(2,1) ./ se));


	# now the Davidson-MacKinnon bootstrap samples using thetahat but substituting in values under H0
	# (x is fixed for the bootstraps, but varies across montecarlo trials)
	# simulate using slope = 0
	theta(2,1) = 0; # set to value under H0
	thetas = zeros(bs_reps,1); # container for bs output
	for i = 1:bs_reps
		conv = 0;
		while conv != 1
			y = ProbitDGP(theta, x);
			data = [y x];
			[thetaboot, junk, conv]  = mle_estimate(theta, data, model, modelargs, control);
			if conv != 1 restarts = restarts + 1; endif # tally up convergence failures

		endwhile
		thetas(i,:) = thetaboot(2,1);
	endfor
	# find percentage of bootstraps greater than abs value of est.
	p2 = sum(abs(thetas) > abs(thetahat(2,1))) / bs_reps;

	results = [restarts slope p1 p2];

endfunction
