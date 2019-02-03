// The two shock model
// This takes steady state of hours as given, and computes the psi
// parameter that corresponds.
close all;

// Define variables
var y c k n invest z1 z2 MUC MUL r w;
varexo e1 e2;

// Parameters
parameters alppha betta delta gam  nss rho1 sigma1 rho2 sigma2 psi c1 iss yss kss css;
load parameterfile;
set_param_value('alppha',0.33)
set_param_value('betta',beta)
set_param_value('delta',0.025)
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
  r = alppha * exp(z1) * (k(-1))^(alppha-1) * n^(1-alppha);
  w = (1-alppha)*exp(z1)* (k(-1))^alppha * n^(-alppha);
  MUC = betta*MUC(+1) * (1 + r(+1) - delta);
  MUL/MUC = w;
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
stoch_simul(nograph, noprint, nomoments, order=3, periods=260, irf=0) y c n r w;
