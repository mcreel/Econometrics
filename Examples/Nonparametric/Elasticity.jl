using Econometrics

# sample size
n = 30

# draw n values of x between 0 and 4
x = 0:4/(n-1):4

# true function (unknown to researcher)
function f(x)
    10 + 3*x - x^2
end

# derivative of the function
function fprime(x)
    3 - 2x
end

# elasticity of the function: marginal function divided by average function
elasticity = x -> fprime(x)/(f(x)/x)
e = elasticity.(x) # elasticities at the x values


# noisy observation of true function
y = f.(x) + randn(n) 
# researcher only observes y and x


# linear fit and elasticity according to linear fit
X = [ones(n) x]
Xprime = [zeros(n) ones(n)]
b, fit, junk = lsfit(y,X)
elinear = Xprime*b./(fit./x)

# Fourier fit and elasticity
#X = [ones(n) x cos.(x) sin.(x) cos.(2x) sin.(2x)]
#Xprime = [zeros(n) ones(n) -sin.(x) cos.(x) -2sin.(2x) 2cos.(2x)]
X = [ones(n) x cos.(x) sin.(x)]
Xprime = [zeros(n) ones(n) -sin.(x) cos.(x)]
b, fourierfit, junk = lsfit(y,X)
efourier = Xprime*b./(fourierfit./x)


# plot function and fits
p1 = plot(x,[f.(x) fit fourierfit], title = "Function and fits", legend=:bottomleft, label=["true" "linear" "fourier"])

# plot true elasticities and fits
p2 = plot(x, [e elinear efourier], title = "Elasticity and fitted elasticities", legend=:bottomleft, label=["true" "linear" "fourier"])


plot(p1, p2)


