# GMM estimation for a sample from Chi^2(theta)
# compare to two method of moments estimators (see chi2mm.m)
using Distributions, Plots
function main()
n = 30
theta = 3.0
reps = 1000
results = zeros(reps)
W = eye(2)
for i = 1:1000
    y = rand(Chisq(theta), n)
    obj = theta -> ([theta.-mean(y); theta.-0.5*var(y)]'W*[theta.-mean(y); theta.-0.5*var(y)]) 
    thetahat, junk, junk = fminunc(obj, [3.0])
    results[i] = sqrt(n)*(thetahat[1]-theta)
end    
histogram(results,nbins=50)
gui()
return
end
main()
