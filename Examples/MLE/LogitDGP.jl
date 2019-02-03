#  This returns data that follow a logit model with the
#  supplied parameters. n is the number of observations.
#  Note how the binary 0/1 variable is generated so that
#  it really follows the logit model
function LogitDGP(n, theta)
	k = size(theta,1)
	x = ones(n,1)
	if k>1 
		x = [x  randn(n,k-1)]
	end
	y = Float64.((1.0 ./ (1.0 .+ exp.(-x*theta)) .> rand(n,1)))
    return y, x
end


