function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
% function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
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
%   g2
%

if T_flag
    T = CKml.dynamic_g2_tt(T, y, x, params, steady_state, it_);
end
v2 = zeros(34,3);
v2(1,1)=1;
v2(2,1)=2;
v2(3,1)=3;
v2(4,1)=3;
v2(5,1)=3;
v2(6,1)=3;
v2(7,1)=3;
v2(8,1)=3;
v2(9,1)=3;
v2(10,1)=3;
v2(11,1)=3;
v2(12,1)=4;
v2(13,1)=4;
v2(14,1)=4;
v2(15,1)=4;
v2(16,1)=4;
v2(17,1)=4;
v2(18,1)=4;
v2(19,1)=4;
v2(20,1)=4;
v2(21,1)=5;
v2(22,1)=5;
v2(23,1)=6;
v2(24,1)=6;
v2(25,1)=6;
v2(26,1)=9;
v2(27,1)=9;
v2(28,1)=9;
v2(29,1)=9;
v2(30,1)=9;
v2(31,1)=9;
v2(32,1)=9;
v2(33,1)=9;
v2(34,1)=9;
v2(1,2)=77;
v2(2,2)=172;
v2(3,2)=1;
v2(4,2)=7;
v2(5,2)=109;
v2(6,2)=9;
v2(7,2)=145;
v2(8,2)=115;
v2(9,2)=117;
v2(10,2)=151;
v2(11,2)=153;
v2(12,2)=1;
v2(13,2)=7;
v2(14,2)=109;
v2(15,2)=9;
v2(16,2)=145;
v2(17,2)=115;
v2(18,2)=117;
v2(19,2)=151;
v2(20,2)=153;
v2(21,2)=268;
v2(22,2)=285;
v2(23,2)=191;
v2(24,2)=192;
v2(25,2)=209;
v2(26,2)=1;
v2(27,2)=7;
v2(28,2)=109;
v2(29,2)=9;
v2(30,2)=145;
v2(31,2)=115;
v2(32,2)=117;
v2(33,2)=151;
v2(34,2)=153;
v2(1,3)=(-(getPowerDeriv(y(5),(-params(4)),2)));
v2(2,3)=(-(params(10)*exp(y(10))));
v2(3,3)=(-(T(2)*params(1)*exp(y(9))*getPowerDeriv(y(1),params(1)-1,2)));
v2(4,3)=(-(T(5)*T(7)));
v2(5,3)=v2(4,3);
v2(6,3)=(-(T(2)*T(5)));
v2(7,3)=v2(6,3);
v2(8,3)=(-(T(1)*T(9)));
v2(9,3)=(-(T(1)*T(7)));
v2(10,3)=v2(9,3);
v2(11,3)=(-(T(1)*T(2)));
v2(12,3)=(-(T(4)*exp(y(9))*(1-params(1))*T(10)));
v2(13,3)=(-(exp(y(9))*(1-params(1))*T(6)*T(8)));
v2(14,3)=v2(13,3);
v2(15,3)=(-(T(4)*exp(y(9))*(1-params(1))*T(6)));
v2(16,3)=v2(15,3);
v2(17,3)=(-(exp(y(9))*(1-params(1))*T(3)*getPowerDeriv(y(7),(-params(1)),2)));
v2(18,3)=(-(exp(y(9))*(1-params(1))*T(3)*T(8)));
v2(19,3)=v2(18,3);
v2(20,3)=(-(exp(y(9))*(1-params(1))*T(3)*T(4)));
v2(21,3)=(-params(2));
v2(22,3)=v2(21,3);
v2(23,3)=(-((-y(12))*(y(11)+y(11))))/(y(11)*y(11)*y(11)*y(11));
v2(24,3)=(-1)/(y(11)*y(11));
v2(25,3)=v2(24,3);
v2(26,3)=(-(T(2)*exp(y(9))*T(10)));
v2(27,3)=(-(exp(y(9))*T(6)*T(7)));
v2(28,3)=v2(27,3);
v2(29,3)=(-(T(2)*exp(y(9))*T(6)));
v2(30,3)=v2(29,3);
v2(31,3)=(-(exp(y(9))*T(3)*T(9)));
v2(32,3)=(-(exp(y(9))*T(3)*T(7)));
v2(33,3)=v2(32,3);
v2(34,3)=(-(T(2)*exp(y(9))*T(3)));
g2 = sparse(v2(:,1),v2(:,2),v2(:,3),11,324);
end
