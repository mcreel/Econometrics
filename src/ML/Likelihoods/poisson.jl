function  poisson(theta, y, x)
	lambda = exp.(x*theta)
    logdensity = y .* (x*theta) .- lambda .- logabsgamma.(y.+1.0)[1]
end
