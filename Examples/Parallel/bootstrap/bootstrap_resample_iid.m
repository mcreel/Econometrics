# take a bootstrap sample from a data matrix
# if second arg is omitted, the sample has same size as original
function newsample = bootstrap_resample_iid(data, reps)
	n = rows(data);
	if (nargin < 2) reps = n; endif
	index = 1 + (n-1)*rand(reps,1);
	index = round(index);
	newsample = data(index,:);
endfunction

# This is a simple check that it works
# x = 1:10;
# x = x';
# for i = 1:20
# bootstrap(x)'
# endfor

