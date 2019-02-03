# outfile = "Cluster_50_nodes-";

# slaves go from 0 to 10, then up to 50 in increments of 5
#nslaves = [0:10, 5*(1:8)+10];

nslaves = 0:4;
nslaves = nslaves';

nslaves = [0:4]';

# this section defines the data
load data;
k = columns(data) - 1;
model = "NegBinSNP";
modelargs = {1}; # this is constant from now on
theta = zeros(k+4,1);

# bfgs controls: maxiters, verbosity, conv, minarg
control = {15,2,1,1};	

results = zeros(rows(nslaves),3);

for i = 1:rows(nslaves);
	tic();
	[thetahat, junk, junk, iters] = mle_estimate(theta, data, model, modelargs, control, nslaves(i), 0);
	t = toc();
	results(i,:) = [nslaves(i), t, iters];
	printf("nodes %d   time %f  iters %f  \n", nslaves(i), t, iters);
endfor

# eval (sprintf("save \"%sML-NBSNP-results.out\" results", outfile));
