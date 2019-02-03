function output = pois_bootstrap(f_arg)
	data = f_arg{1};
	theta = f_arg{2};
	model = f_arg{3};
	theta = mle_estimate(theta, data, model, {}, {100,0});
	output = theta'; # we want a row vector!
endfunction	
	 
