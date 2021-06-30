function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
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
%   residual
%

if T_flag
    T = CKmcmc.static_resid_tt(T, y, x, params);
end
residual = zeros(11, 1);
lhs = y(8);
rhs = y(2)^(-params(4));
residual(1) = lhs - rhs;
lhs = y(9);
rhs = params(10)*exp(y(7));
residual(2) = lhs - rhs;
lhs = y(10);
rhs = T(1)*T(2);
residual(3) = lhs - rhs;
lhs = y(11);
rhs = exp(y(6))*(1-params(1))*T(3)*T(4);
residual(4) = lhs - rhs;
lhs = y(8);
rhs = y(8)*params(2)*(1+y(10)-params(3));
residual(5) = lhs - rhs;
lhs = y(9)/y(8);
rhs = y(11);
residual(6) = lhs - rhs;
lhs = y(6);
rhs = y(6)*params(6)+params(7)*x(1);
residual(7) = lhs - rhs;
lhs = y(7);
rhs = y(7)*params(8)+params(9)*x(2);
residual(8) = lhs - rhs;
lhs = y(1);
rhs = T(2)*exp(y(6))*T(3);
residual(9) = lhs - rhs;
lhs = y(5);
rhs = y(1)-y(2);
residual(10) = lhs - rhs;
lhs = y(3);
rhs = y(5)+y(3)*(1-params(3));
residual(11) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
