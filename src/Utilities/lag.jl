# returns the variable (or matrix), lagged p times,
# with the first p rows filled with ones (to avoid divide errors)
# remember to drop those rows before doing analysis
function lag(x::Array{Float64,2},p::Int64)
	n,k = size(x)
	lagged_x = [ones(p,k); x[1:n-p,:]]
end

function lag(x::Array{Float64,1},p::Int64)
	n = size(x,1)
	lagged_x = [ones(p); x[1:n-p]]
end	 
