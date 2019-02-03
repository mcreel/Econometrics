
%{
Set the true parameters within these bounds
   Lower     Upper 
   0.20000   0.40000
   0.95000   1.00000
   0.01000   0.10000
   0.00000   5.00000
   0.00000   1.00000
   0.00000   0.10000
   0.00000   1.00000
   0.00000   0.10000
   0.25000   0.37500
%}
% set true parameters here
alpha  = 0.3;
beta   = 0.99;
delta  = 0.025;
gam    = 3.5;
rho1   = 0.5;
sigma1 = 0.03;
rho2   = 0.3;
sigma2 = 0.01;
nss    = 1/3;
theta0 = [alpha; beta; delta; gam; rho1; sigma1; rho2; sigma2; nss];

% the psi implied by other parameters
c1 = ((1/beta + delta - 1)/alpha)^(1/(1-alpha));
kss = nss/c1;
css = kss * (c1^(1-alpha) - delta);
c2 = (css)^(-gam/alpha);
psi = (1-alpha)*((c1/c2)^(-alpha));

save parameterfile  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
RNGstate = rand('state');   % this is used to get different random draws
                            % each time, otherwise, Dynare fixes the state
                            % so the results are the same each time this is run
dynare SimpleModel noclearall; % solve the model
% set a new state in case this is run again
ss = round(RNGstate(1)/RNGstate(2)*sum(RNGstate));
set_dynare_seed(ss);
ok = false;
while !ok
    % do the simulation at the param values
    info = stoch_simul(var_list_);
    % get a simulation of length 160 and compute aux. statistic
    data = [y c n MPK MPL];
    data = data(101:260,:);
    Z = aux_stat(data);
    ok = Z(1,1) != -1000;
endwhile
% compute the output of the net at the input values of the statistics
thetahat = NNstat(Z');
fprintf("   true       estimated\n");
disp([theta0 thetahat']);
plot(data);
system("./cleanup");
