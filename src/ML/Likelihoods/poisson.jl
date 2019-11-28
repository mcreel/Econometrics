function  poisson(theta, y, x)
	lambda = exp.(x*theta)
    n = size(y,1)
    c = zeros(n)
    for i = 1:n
        c[i] = logabsgamma(y[i])[1]
    end    
    logdensity = y .* (x*theta) .- lambda .- c
end
