function [data, L, K, snm_draws] = ar1_model(theta, modelargs)
	snm_draws = modelargs{1};
	simlength = modelargs{2};
	burnin = modelargs{3};
	make_instruments = modelargs{4};
	if ischar(snm_draws)
		snm_draws = randn(simlength + burnin + 10, 2);
		modelargs{1}= snm_draws;
	endif

	n = modelargs{2};
	maxlag = 1;
	n = n + maxlag;
	modelargs{2} = n;

	y = ar1(theta, modelargs);
	data = [y lags(y,maxlag)];
	data = data(maxlag + 1:n,:);
	n = rows(data);

	y = data(:,1);
	endogs = y;
	condvars = data(:,2);
	data = [endogs condvars];
	L = columns(endogs);
	K = columns(condvars);

	if make_instruments
		instruments = [ones(n,1) st_norm(condvars)];
		data = [data instruments];
	endif

endfunction
