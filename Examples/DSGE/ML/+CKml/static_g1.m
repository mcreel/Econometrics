function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = CKml.static_g1_tt(T, y, x, params);
end
g1 = zeros(11, 11);
g1(1,2)=(-(getPowerDeriv(y(2),(-params(4)),1)));
g1(1,8)=1;
g1(2,7)=(-(params(10)*exp(y(7))));
g1(2,9)=1;
g1(3,3)=(-(T(2)*params(1)*exp(y(6))*getPowerDeriv(y(3),params(1)-1,1)));
g1(3,4)=(-(T(1)*T(6)));
g1(3,6)=(-(T(1)*T(2)));
g1(3,10)=1;
g1(4,3)=(-(T(4)*exp(y(6))*(1-params(1))*T(5)));
g1(4,4)=(-(exp(y(6))*(1-params(1))*T(3)*getPowerDeriv(y(4),(-params(1)),1)));
g1(4,6)=(-(exp(y(6))*(1-params(1))*T(3)*T(4)));
g1(4,11)=1;
g1(5,8)=1-params(2)*(1+y(10)-params(3));
g1(5,10)=(-(y(8)*params(2)));
g1(6,8)=(-y(9))/(y(8)*y(8));
g1(6,9)=1/y(8);
g1(6,11)=(-1);
g1(7,6)=1-params(6);
g1(8,7)=1-params(8);
g1(9,1)=1;
g1(9,3)=(-(T(2)*exp(y(6))*T(5)));
g1(9,4)=(-(exp(y(6))*T(3)*T(6)));
g1(9,6)=(-(T(2)*exp(y(6))*T(3)));
g1(10,1)=(-1);
g1(10,2)=1;
g1(10,5)=1;
g1(11,3)=1-(1-params(3));
g1(11,5)=(-1);
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
