# overall controls for experiment
n = 15; # sample size
iter_limit = 100; # limit BFGS iterations, if hit, the draw is not used
n_pooled = 1; # slaves do this many reps before sending results back to master

# SIZE
slope = 0; # true slope
outfile = "size.out"; # where to write results
bs_reps = 399; # number of boostrap replications

	# half on big cluster
	reps = 2500; # number of monte carlo replications
	nslaves = 19; # how many slave nodes?
	montecarlo("bs_test", {n, bs_reps, slope, iter_limit}, reps, outfile, nslaves, n_pooled);

	# half on small cluster
	reps = 2500; # number of monte carlo replications
	nslaves = 9; # how many slave nodes?
	montecarlo("bs_test", {n, bs_reps, slope, iter_limit}, reps, outfile, nslaves, n_pooled);

# POWER
slope = 1; # true slope
outfile = "power.out"; # where to write results
bs_reps = 999; # number of boostrap replications

	# half on big cluster
	reps = 1000; # number of monte carlo replications
	nslaves = 19; # how many slave nodes?
	montecarlo("bs_test", {n, bs_reps, slope, iter_limit}, reps, outfile, nslaves, n_pooled);

	# half on small cluster
	reps = 1000; # number of monte carlo replications
	nslaves = 9; # how many slave nodes?
	montecarlo("bs_test", {n, bs_reps, slope, iter_limit}, reps, outfile, nslaves, n_pooled);
