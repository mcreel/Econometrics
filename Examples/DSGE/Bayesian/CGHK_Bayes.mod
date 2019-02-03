// The two shock model of the paper 'Indirect Likelihoood Inference'
// by Michael Creel and Dennis Kristensen
// This takes steady state of hours as given, and computes the psi
// parameter that corresponds.
close all;

// Define variables
var y c k n invest z1 z2 MUC MUL MPK MPL;
varexo e1 e2;

// Parameters
parameters alppha betta delta gam  nss rho1 sigma1 rho2 sigma2 psi c1 iss yss kss css;
load parameterfile;
set_param_value('alppha',alpha)
set_param_value('betta',beta)
set_param_value('delta',delta)
set_param_value('gam',gam)
set_param_value('rho1',rho1)
set_param_value('sigma1',sigma1)
set_param_value('rho2',rho2)
set_param_value('sigma2',sigma2)
set_param_value('nss',nss)

// Model
model;
  MUC = (c)^(-gam);
  MUL = psi*exp(z2);
  MPK = alppha * exp(z1) * (k(-1))^(alppha-1) * n^(1-alppha);
  MPL = (1-alppha)*exp(z1)* (k(-1))^alppha * n^(-alppha);
  MUC = betta*MUC(+1) * (1 + MPK(+1) - delta);
  MUL/MUC = MPL;
  z1 = rho1*z1(-1) + sigma1*e1;
  z2 = rho2*z2(-1) + sigma2*e2;
  y = exp(z1) * ((k(-1))^alppha) * (n^(1-alppha));
  invest = y - c;
  k = invest + (1-delta)*k(-1);
end;


// Shocks
shocks;
var e1 = 1;
var e2 = 1;
end;

% note: estimating alpha, gam, rho1, sig1, rho2, sig2 and nss works
% if obs. vars. are c and n
% trying to estimate other parameters (e.g., beta, delta) shows problems

estimated_params;
betta, uniform_pdf, , 0.99,  0.95, 0.9999; 
gam, uniform_pdf, , 2.0, 0.0,  5.0;
rho1, uniform_pdf, , 0.9, 0.0,  1.0; 
sigma1, uniform_pdf, , 0.02, 0.0,  0.1;
rho2, uniform_pdf, , 0.7,   0.0,  1.0; 
sigma2, uniform_pdf, , 0.01, 0.0,  0.1;
nss, uniform_pdf, , 0.33, 6/24, 9/24; 
end;

varobs c n;  // experiment choosing one or two from y c n MPK MPL (c and n work well)

// short MCMC for illustration only 
// set order=2 to see particle filter option
estimation(datafile=dsgedata, nobs=160, order=1, mh_replic=10000, mh_nblocks=1, mh_jscale=0.8, mh_init_scale=5, irf=40) y c n MPK MPL;



