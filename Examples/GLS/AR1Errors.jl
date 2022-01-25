# this script does a Monte Carlo that shows efficiency of GLS compared to OLS
# when there are AR1 errors
using Econometrics, Plots

# the function to be Monte Carlo'ed
function wrapper(args)
	n, ρ = args    # sample size and AR parameter
	x = randn(n)     # an exogenous regressor
	# make the AR1 errors
    ϵ = zeros(n)
    ϵ[1] = randn()/sqrt(1-ρ^2)
	for t = 2:n
        ϵ[t] = ρ*ϵ[t-1] + randn()
    end
	# the dep var
	y = 1. .+ x .+ ϵ
	# OLS
	x = [ones(n) x]
	β_ols = x \ y
	# estimate ρ
	e = y - x*β_ols
    # estimating ρ using correlation, to keep in stationary region
    ρhat = cor(e[1:end-1], e[2:end]) 
	ystar = y - ρhat.*lag(y,1)
	xstar = x - ρhat*lag(x,1)
    ystar[1] = sqrt(1-ρhat^2) .*y[1]
    xstar[1,:] = sqrt(1-ρhat^2) .* x[1,:]
	β_gls = xstar\ystar
	β_ols, β_gls
end

# do the Monte Carlo
n = 30
ρ = 0.9
args = (n, ρ)
reps = 1000
β_ols = zeros(reps,2)
β_gls = zeros(reps,2)
for i = 1:reps
    β_ols[i,:], β_gls[i,:] = wrapper(args)
end    

# analyze results
p1 = histogram(β_ols[:,2],bins=30, legend=false, xlabel="OLS")
p2 = histogram(β_gls[:,2],bins=30, legend=false, xlabel="GLS")
plot(p1,p2)
#savefig("AR1errors_OLSvsGLS.png")
gui()
		

