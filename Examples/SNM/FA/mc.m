if !(exist('./fa.oct','file')) system("mkoctfile fa.cc"); endif
n = 500; # number of obs
outfile = "fa.out";
model = "fa_model"; # the DGP
# true parameters
a1 = 0.2;
a2 = 0.7;
sig = 0.5;
b2 = -0.5;
theta = [a1; a2; sig; b2];

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

