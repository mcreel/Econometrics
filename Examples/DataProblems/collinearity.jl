# illustrates effects of collinearity
# try setting rho=0.1 or rho=0.9, and observe
# the effect on standard dev of coefficients
using Econometrics, StatsPlots, LinearAlgebra
function wrapper(rho, n)
  V = [1.0 rho; rho 1.0]
  P = cholesky(V).U
  x = randn(n,2)*P
  x = [ones(n,1) x]
  b = ones(3)
  y = x*b + randn(n,1)
  b1 = x\y

  z = x[:,1:2]
  b2 = z\y

  R = [0.0 1.0 -1.0]
  r = 0.0
  b3 = ols(y, x, R=R, r=r, silent=true)[1]
  results = [b1' b2' b3']
end

function main()
reps = 1000
rho = 0.9
n = 30
data = zeros(reps,8)
for rep = 1:reps
    data[rep,:] = wrapper(rho,n)
end
println("correlation between x2 and x3: ",rho)
println("descriptive statistics for 1000 OLS replications")
dstats(data[:,1:3])
println("descriptive statistics for 1000 OLS replications, dropping x3")
dstats(data[:,4:5])
println("descriptive statistics for 1000 Restricted OLS replications, b2=b3")
dstats(data[:,6:8])
density(data[:,2],label="OLS")
density!(data[:,5],label="Drop X3")
density!(data[:,7],label="Restricted LS", title="collinearity: OLS and restricted OLS - true β₂=1")
gui()
savefig("collin.png")
return
end
main()
