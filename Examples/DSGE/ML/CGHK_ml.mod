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

// generate data
%stoch_simul(nograph, noprint, nomoments, order=3, periods=260, irf=0) y c n MPK MPL;

% note: estimating alpha, beta, gam,  rho1, sig1, rho2, sig2 and nss works
% if obs. vars. are c and n
% trying to estimate other parameters (delta) shows problems

estimated_params;
// the numbers after the comment are the true 
// values used to generate the data. Try setting nss=0.2, gam=1, and rho1=0.5
betta, 0.99;    // 0.99
gam, 2.0;       // 2.0
rho1, 0.9;      // 0.9
sigma1, 0.02;   // 0.02
rho2, 0.7;      // 0.7
sigma2, 0.01;   // 0.01
nss, .333;       // 1/3
end;

# using c and n gives good results, c and MPL, not so good
varobs n c;  // experiment choosing one or two from y c n MPK MPL

// computes only the posterior mode for demonstration. 
estimation(datafile=dsgedata, nobs=160, order=1);

