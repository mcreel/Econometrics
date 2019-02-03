# number of monte carlo replications
reps = 1000;

# specifics of test
T = 1000;
dim = 5;

# run it on master and slave(s)
compute_nodes = 0;
n_pooled = 10;
debug = false;
close;
montecarlo("tracetest", {T, dim}, reps, "example1.out", compute_nodes, n_pooled);

# examine results
load example1.out; # load output
output = example1(:,3); # get results of tracetest
hist(output, 30);
xlabel("");
title("Monte Carlo example, parallel");
