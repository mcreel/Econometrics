function m = Poisson_IV_moments(theta, data, momentargs)
	k = momentargs{1}; # where does x end and w starts?
	y = data(:,1);
	x = data(:,2:k+1);
	w = data(:, k+2:columns(data));
	lambda = exp(x*theta);
	e = y  ./ lambda - 1;
	m = dmult(e, w);
endfunction	 
