# we need the input argument to be a cell, so this is 
# a little wrapper function for ols()
function result = myols(f_args)
	data = f_args{1};
	k = columns(data) - 1;
	y = data(:,1);
	x = data(:,2:k+1);
	beta = ols(y,x);
	result = beta';
endfunction	
