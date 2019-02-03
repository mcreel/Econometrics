# ridge regression as example of Bayesian estimation
using Plots, LinearAlgebra
function main()
n = 100
a = 0.95 # controls collinearity: close to 0 indep close to 1 highly collinear
bs = zeros(1000,6)
truebeta = ones(3,1)
for i = 1:1000
	# generate regressors that are collinear
	x1 = rand(n,1)
	x2 = a*x1 + (1-a)*randn(n,1)
	x = [ones(n,1) x1 x2]
	# true params and dep var
	e = randn(n,1)
	y = x*truebeta+e;
	# OLS
	bols = x\y
	# ridge regression
	k = 1.0
	yy = [y; k*zeros(3,1)]
	xx = [x; k*Matrix{Float64}(I, 3, 3)]
	bridge = xx\yy
	# store results
	bs[i,:] = [bols' bridge']

end
# histograms to compare 
bs = bs .- [truebeta' truebeta'] # take away true, so centered around 0
p1 = histogram(bs[:,2], title="OLS")
p2 = histogram(bs[:,5],title="Ridge")
plot(p1,p2,layout=(1,2), legend=false)
#savefig("RidgeExample.svg")
gui()
return
end
main()
