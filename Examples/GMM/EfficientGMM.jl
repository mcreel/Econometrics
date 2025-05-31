# GMM estimation for a sample from Chi^2(theta)
# compare to two method of moments estimators (see chi2mm.m)
using Econometrics, Random, Distributions, Plots
function main()
n = 30
theta = 3.0
reps = 1000
results = zeros(reps,2)

y = zeros(n) # this is just a place holder to define the moments, below
# define the moment conditions
m = theta -> [(theta.- y) (theta .- 0.5*(y .- theta).^2.0)] 
for i = 1:1000
    rand!(Chisq(theta), y) # this replaces the place holder, from above
    thetahat, junk, junk, ms, junk = gmm(m, [3.0], eye(2))
    results[i,1] = sqrt(n)*(thetahat[1]-theta)
    W = inv(cov(ms))
    thetahat, junk, ms, junk = gmm(m, thetahat, W)
    results[i,2] = sqrt(n)*(thetahat[1]-theta)
end    
p1 = npdensity(results[:,1])
plot!(p1, title="Inefficient")
p2 = npdensity(results[:,2])
plot!(p2, title="Efficient")
p = plot(p1,p2,layout=(2,1),xlims=(-10,20))
#savefig("Efficient.png")
println("Monte Carlo covariance: (1,1) is inefficient GMM, (2,2) is efficient")
prettyprint(cov(results))
return p
end
main()

