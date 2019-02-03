% returns the variable (or matrix), lagged p times,
% with the first p rows filled with ones (to avoid divide errors)
% remember to drop those rows before doing analysis
function lagged_x = lag(x,p)
	k = size(x,2);
	lagged_x = [ones(p,k); x(1:size(x,1)-p,:)];
end	 
