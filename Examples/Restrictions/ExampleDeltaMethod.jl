# This file illustrates how the delta method may be used
# to calculate the covariance of a nonlinear function of
# estimated parameters
using ForwardDiff LinearAlgebra, Statistics
# this function that gives the elasticities of x*b wrt x:
function ElasticitiesLinearModel(theta, x)
	elasticities = (theta .* x) / (x'theta)
end

function main()
# define some random data	
n = 100
k = 5
x = [ones(n,1) rand(n,k-1)] 
beta = (1:k) # true betas are [1,2,3,4,5]
y = x*beta + 3.0*randn(n)
b, varb, junk = ols(y,x)

# we'll evaluate elasticities at the sample means of the data
x = vec(mean(x,dims=1))

Elasticities = theta -> vec(ElasticitiesLinearModel(theta, x))
# want derivatives by param in the rows: kXg from GMM theory 
J = ForwardDiff.jacobian(Elasticities, b)'
elasticities = Elasticities(b)
var_elasticities = J*varb*J'
# Taa daa - the results!
println()
println("Delta method elasticities and standard errors, at the means of the regressors")
prettyprint([elasticities diag(var_elasticities)])
return
end
main()
