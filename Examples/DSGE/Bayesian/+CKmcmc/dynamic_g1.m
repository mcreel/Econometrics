function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = CKmcmc.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(11, 18);
g1(1,5)=(-(getPowerDeriv(y(5),(-params(4)),1)));
g1(1,11)=1;
g1(2,10)=(-(params(10)*exp(y(10))));
g1(2,12)=1;
g1(3,1)=(-(T(2)*params(1)*exp(y(9))*getPowerDeriv(y(1),params(1)-1,1)));
g1(3,7)=(-(T(1)*T(6)));
g1(3,9)=(-(T(1)*T(2)));
g1(3,13)=1;
g1(4,1)=(-(T(4)*exp(y(9))*(1-params(1))*T(5)));
g1(4,7)=(-(exp(y(9))*(1-params(1))*T(3)*getPowerDeriv(y(7),(-params(1)),1)));
g1(4,9)=(-(exp(y(9))*(1-params(1))*T(3)*T(4)));
g1(4,14)=1;
g1(5,11)=1;
g1(5,15)=(-(params(2)*(1+y(16)-params(3))));
g1(5,16)=(-(params(2)*y(15)));
g1(6,11)=(-y(12))/(y(11)*y(11));
g1(6,12)=1/y(11);
g1(6,14)=(-1);
g1(7,2)=(-params(6));
g1(7,9)=1;
g1(7,17)=(-params(7));
g1(8,3)=(-params(8));
g1(8,10)=1;
g1(8,18)=(-params(9));
g1(9,4)=1;
g1(9,1)=(-(T(2)*exp(y(9))*T(5)));
g1(9,7)=(-(exp(y(9))*T(3)*T(6)));
g1(9,9)=(-(T(2)*exp(y(9))*T(3)));
g1(10,4)=(-1);
g1(10,5)=1;
g1(10,8)=1;
g1(11,1)=(-(1-params(3)));
g1(11,6)=1;
g1(11,8)=(-1);

end
