# GMM estimation for a sample from Chi^2(theta)
# compare to two method of moments estimators (see chi2mm.m)
using Econometrics, Random, LinearAlgebra, Distributions, Plots
function main()
n = 30
theta = 3.0
reps = 1000
results = zeros(reps)
W = eye(2)
y = zeros(n) # this is just a place holder to define the moments, below
# define the moment conditions
m = theta -> [theta.-mean(y); theta.-0.5*var(y)] 
for i = 1:1000
    rand!(Chisq(theta), y) # this replaces the place holder, from above
    obj = theta -> dot(m(theta),W*m(theta)) 
    thetahat, junk, junk = fminunc(obj, [3.0])
    results[i] = sqrt(n)*(thetahat[1]-theta)
end    
p = histogram(results,nbins=50)
end
main()
