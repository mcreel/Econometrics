# GMM estimation for a sample from Chi^2(theta)
# compare to two method of moments estimators (see chi2mm.m)
using Distributions, Plots
function main()
n = 30
theta = 3.0
reps = 1000
results = zeros(reps,2)
for i = 1:1000
    y = rand(Chisq(theta), n)
    moments = theta -> [theta.-y  theta.-0.5.*(y .- mean(y)).^2.0]
    thetahat, junk, junk, ms, junk = gmm(moments, [3.0], eye(2))
    results[i,1] = sqrt(n)*(thetahat[1]-theta)
    W = inv(cov(ms))
    thetahat, junk, ms, junk = gmm(moments, thetahat, W)
    results[i,2] = sqrt(n)*(thetahat[1]-theta)
end    
p1 = npdensity(results[:,1])
plot!(p1, title="Inefficient")
p2 = npdensity(results[:,2])
plot!(p2, title="Efficient")
plot(p1,p2,layout=(2,1),xlims=(-10,20))
savefig("Efficient.svg")
gui()
println("Monte Carlo covariance: (1,1) is inefficient GMM, (2,2) is efficient")

prettyprint(cov(results))
return
end
main()

