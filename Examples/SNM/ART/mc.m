if !(exist('./art.oct','file')) system("mkoctfile art.cc"); endif
n = 150; # number of obs
outfile = "art.out";
model = "art_model"; # the DGP
# true parameters
a = 0;
b = 0.5;
sig = 1;
theta = [a; b; sig];


S = 10000;  # number of simulations
burnin = 100;

n_pooled = 1;
vc_reps = 100;

#----------------  now sample at true value  ---------------------#
reps = 500;
from_param_space = false;
wrapper_args = {model, theta, n, S, burnin, vc_reps, from_param_space};
#if not(MPI_Initialized) MPI_Init; endif
montecarlo("wrapper", wrapper_args, reps, outfile,  n_pooled, true);
#if not(MPI_Finalized) MPI_Finalize; endif

