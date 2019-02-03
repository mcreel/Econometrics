function  poisson(theta, y, x)
	lambda = exp.(x*theta)
	logdensity = y .* (x*theta) .- lambda .- lgamma.(y.+1.0)
end
