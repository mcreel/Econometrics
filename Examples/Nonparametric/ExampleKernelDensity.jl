n = 30000 # sample size
x = sum(randn(n,3).^2,dims=2) # generates 10000 draws from χ²(3)
using StatsPlots, Distributions
density(x, label="nonparametric density") # creates a nonparametric density plot from the sample
plot!(Chisq(3),label="χ²(3) density") # plots the true density on top
