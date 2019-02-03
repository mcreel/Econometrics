% returns the variable (or matrix), lagged from 1 to p times,
% with the first p rows filled with ones (to avoid divide errors)
% remember to drop those rows before doing analysis
function lagged_x = lags(x,p)
	k = size(x,2);
	lagged_x = lag(x,1);
	if p > 1
		for i = 2:p
			lagged_x = [lagged_x, lag(x,i)]; 
		end
	end	
end	 
