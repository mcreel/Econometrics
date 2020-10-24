# this shows effects of het. on the OLS estimator
using Plots, Distributions, Econometrics, LinearAlgebra, Statistics
function main()
het = true  # set this to true or false
n = 50 # sample size

# the function we simulate
function wrapper(het, n)
  x = randn(n)
  e = randn(n)
  if het
    e = x .* e # note the heteroscedasticity here
  end
  y = x + e # true coefs are zero for const, and 1 for slope
  x = [ones(n) x]
  b, varb, junk, junk, junk = ols(y, x, vc="standard",silent=true)
  H0 = [0; 1] # coefficient values under true null hypothesis
  t = (b - H0) ./ sqrt.(diag(varb))
end

# simulation loop
reps = 1000
ts = zeros(reps,2)
for rep = 1:reps
    ts[rep,:] = wrapper(het, n)
end    
ts = abs.(ts)
crit_val = quantile(TDist(n-2),0.95) # 5# of prob. is to the R of this value
test = ts .> crit_val # now it"s a 10# signif. level test, because of the abs value
if het
  println("data is heteroscedastic")
else
  println("data is homoscedastic")
end
println("rejection frequency of nominal 10 percent test")
println("intercept: ", sum(test[:,1]/reps))
println("slope: ", sum(test[:,2]/reps))

plot([0.0; 1.0],sum(test,dims=1)'/1000, linetype=:bar, label="0 is intercept, 1 is slope")
#savefig("EffectsOLS.png")
gui()
return
end
main()
