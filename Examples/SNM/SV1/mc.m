if !(exist('./sv.oct','file')) system("mkoctfile sv.cc"); endif
n = 500; # number of obs
outfile = "sv.out";
model = "sv_model"; # the DGP

# true parameters
phi = 0.9;
sig_b = exp(-0.736/2);
sig_e = 0.363;
theta = [phi; sig_b; sig_e];

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

