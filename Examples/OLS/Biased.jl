# this is a little Monte Carlo exercise that illustrates that
# the OLS estimator is biased when we have an autoregressive model
# with weak exogeneity
using Plots
function main()
reps = 1000 # number of Monte Carlo reps.
n = 20 # sample size
betas = zeros(reps)
truebetas = [0, 0.9]
for i = 1:reps
	x = zeros(n+1)
	x[1] = 0.0
	# generate AR(1) data
	for t = 2:n+1
		x[t] = truebetas[1] + truebetas[2]*x[t-1] + randn()
	end
	y = x[2:n+1]    # dependent variable
	x = x[1:n]      # explanatory variable is the lagged dep var.
	x = [ones(n,1) x]
	beta = y\x
	betas[i] = beta[2]
end
betas = betas .- truebetas[2]
histogram(betas, label="", show=true)
#savefig("Biased.png")
return
end
main()

