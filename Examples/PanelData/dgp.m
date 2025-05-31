function [y, x] = dgp(n, T, k)
	x = randn(n*T,k);  	% regressors
	e = randn(n*T,1);	% error (white noise)
	% agent effect: even grid over (0,2) 
	temp = linspace(0,2,n)';
	agent_effect = kron(temp,ones(T,1));
	% time effect: even grid over (-1,1)
	%temp = linspace(0,2,T)';
	%time_effect = kron(ones(n,1),temp);
	time_effect = 0;
	beta = ones(k,1);	% true slopes
	y = -1 + agent_effect + time_effect + x*beta + e;
end

