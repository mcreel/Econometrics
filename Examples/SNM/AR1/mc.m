if !(exist('./ar1.oct','file')) system("mkoctfile ar1.cc"); endif
n = 50; # number of obs
outfile = "ar1.out";
model = "ar1_model"; # the DGP
# true parameters
theta = [0; 0.9];

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

