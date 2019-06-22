# returns the variable (or matrix), lagged from 1 to p times,
# with the first p rows filled with ones (to avoid divide errors)
# remember to drop those rows before doing analysis
function  lags(x::Array{Float64,2},p)
	n, k = size(x)
	lagged_x = zeros(eltype(x),n,p*k)
	for i = 1:p
		lagged_x[:,i*k-k+1:i*k] = lag(x,i)
	end
    return lagged_x
end	

function  lags(x::Array{Float64,1},p)
	n = size(x,1)
	lagged_x = zeros(eltype(x), n,p)
	for i = 1:p
		lagged_x[:,i] = lag(x,i)
	end
    return lagged_x
end	 
