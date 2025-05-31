n = 20;
burnin = 100;
rho = 0.9;
sig = 1;
bootstrap_reps = 100;

function [rhohat, uhat] = make_rhohat(rho, u, n, burnin) 
	y = zeros(n+burnin,1);
	y(1,:) = 0;
	for t = 2:n+burnin
		y(t,:) = rho*y(t-1,:) + u(t,:);
	endfor
	y = y(burnin+1:end,:);
	x = y(1:end-1,:);
	y = y(2:end,:);
	rhohat = ols(y(2:end,:), y(1:end-1,:));
	if rhohat > 0.999 rhohat=0.999; endif
	uhat = y-rhohat*x;	
endfunction

% real rhohat
u = sig*randn(n+burnin,1);
[rhohat, uhat] = make_rhohat(rho, u, n, burnin);

% bootstrap resamples to get distribution of rhohat
rhob = zeros(bootstrap_reps,1);
for i = 1:bootstrap_reps
	ub = bootstrap_resample_iid(uhat, n + burnin);
	rhob(i,:) = make_rhohat(rhohat, ub, n, burnin);
endfor
rhohat
mean(rhob)


