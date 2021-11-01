using SpecialFunctions
function  poisson(theta, y, x)
    lambda = exp.(x*theta)
    n = size(y,1)
    logdensity = y .* (x*theta) .- lambda .- loggamma.(y+1)
end
