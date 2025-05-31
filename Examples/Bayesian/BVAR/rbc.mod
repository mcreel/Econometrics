% this is Fernández-Villaverde's basic model rbc.mod from the dynare site,
% but with psi computed to make nss = 1/3
% the other parameters are free
% I use n for labor instead if l, which looks like a 1
% also, the variable y_l (y_n) is dropped, since not used

% Basic RBC Model 
%
% Jesus Fernandez-Villaverde
% Philadelphia, March 3, 2005

%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------

close all;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var y c k i n z;
varexo e;

// Parameters
parameters alpha beta delta psi rho sigma nss;
alpha   = 0.33;
beta    = 0.99;
delta   = 0.023;
rho     = 0.95;  
sigma   = (0.007/(1-alpha));
nss = 1/3;


% this is the part that finds psi to make lss equal specified value
phi = ((1/alpha)*(1/beta -1 + delta)) ^ (1/(1-alpha));
omega = phi^(1-alpha) - delta;
kss = nss/phi;
css = omega*kss;
yss = kss^alpha * nss^(1-alpha);
psi = (1 - nss)/css * (1-alpha)*kss^alpha * nss^(-alpha);



%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model(use_dll); 
  (1/c) = beta*(1/c(+1))*(1+alpha*(k^(alpha-1))*(exp(z(+1))*n(+1))^(1-alpha)-delta);
  psi*c/(1-n) = (1-alpha)*(k(-1)^alpha)*(exp(z)^(1-alpha))*(n^(-alpha));
  c+i = y;
  y = (k(-1)^alpha)*(exp(z)*n)^(1-alpha);
  i = k-(1-delta)*k(-1);
  z = rho*z(-1)+e;
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

initval;
  k = kss;
  c = css;
  n = nss;
  y = yss;
  z = 0; 
  e = 0;
end;

shocks;
var e = sigma^2;
end;

steady;

stoch_simul(order = 3, nograph, noprint, drop=100, periods=400) c i n;

