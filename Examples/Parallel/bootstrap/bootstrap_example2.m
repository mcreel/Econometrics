# outfile = "Cluster_50_nodes-";

reps = 100;

# slaves go from 0 to 10, then up to 50 in increments of 5
nslaves = [0:10, 5*(1:8)+10];
nslaves = nslaves';

nslaves = [0:4]';
nslaves = 0;
load poisson_data;
model = "Poisson"; 	# name of function that calculates loglikelihood
theta = zeros(columns(poisson_data)-1,1); # starting values
f = "pois_bootstrap";
f_arg = {poisson_data, theta, model};


results = zeros(rows(nslaves),2);

for i = 1:rows(nslaves);
	tic();
	bootstrap(f, f_arg, reps, nslaves(i));
	t = toc();
	results(i,:) = [nslaves(i), t];
	printf("slaves %d   time %f  \n", nslaves(i), t);
endfor

# eval (sprintf("save \"%s-results.out\" results", outfile));	
	
