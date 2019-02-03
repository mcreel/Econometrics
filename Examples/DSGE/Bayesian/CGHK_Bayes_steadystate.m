function [ys,check] = junk_steadystate(ys,exe)
global M_ lgy_

if isfield(M_,'param_nbr') == 1
NumberOfParameters = M_.param_nbr;
for i = 1:NumberOfParameters
  paramname = deblank(M_.param_names(i,:));
  eval([ paramname ' = M_.params(' int2str(i) ');']);
end
check = 0;
end


%% Enter model equations here
c1 = ((1/betta + delta - 1)/alppha)^(1/(1-alppha));
kss = nss/c1;
iss = delta*kss;
yss = kss^alppha * nss^(1-alppha);
css = yss - iss;
psi =  (css^(-gam)) * (1-alppha) * (kss^alppha) * (nss^(-alppha));

k = kss;
n = nss;
c = css;
y = yss;
invest = iss;
MUC = c^(-gam);
MPK = alppha  * k^(alppha-1) * n^(1-alppha);
MPL = (1-alppha)* (k)^alppha * n^(-alppha);
MUL = MPL*MUC;
z1 = 0;
z2 = 0;
%% end own model equations


for iter = 1:length(M_.params)
  eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end

if isfield(M_,'param_nbr') == 1

if isfield(M_,'orig_endo_nbr') == 1
NumberOfEndogenousVariables = M_.orig_endo_nbr;
else
NumberOfEndogenousVariables = M_.endo_nbr;
end
ys = zeros(NumberOfEndogenousVariables,1);
for i = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(i,:));
  eval(['ys(' int2str(i) ') = ' varname ';']);
end
else
ys=zeros(length(lgy_),1);
for i = 1:length(lgy_)
    ys(i) = eval(lgy_(i,:));
end
check = 0;
end
