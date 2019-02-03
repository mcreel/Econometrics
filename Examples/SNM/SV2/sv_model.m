function [data, L, K, snm_draws] = sv_model(theta, modelargs)
	snm_draws = modelargs{1};
	simlength = modelargs{2};
	burnin = modelargs{3};
	make_instruments = modelargs{4};
	if ischar(snm_draws)
		snm_draws = randn(simlength + burnin + 10, 2);
		modelargs{1}= snm_draws;
	endif

	theta = parameterize(theta);
	n = modelargs{2};
	maxlag = 4;
	n = n + maxlag;
	modelargs{2} = n;

	y = sv(theta, modelargs);
	y = 25*abs(y);
#dstats(y);
#pause
	data = [y lags(y,maxlag)];
	data = data(maxlag + 1:n,:);
	n = rows(data);

	y = data(:,1);
	endogs = [y y.^2 cos(y) sin(y) cos(2*y) sin(2*y) cos(3*y) sin(3*y) cos(4*y) sin(4*y)];
	condvars = [sum(data(:,2:columns(data)),2)];		
	data = [endogs condvars];
	L = columns(endogs);
	K = columns(condvars);

	if make_instruments
		instruments = [ones(n,1) st_norm(condvars)];
		data = [data instruments];
	endif

	# define data and information needed to split it into pieces
endfunction
