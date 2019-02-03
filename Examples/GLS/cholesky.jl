using Econometrics, LinearAlgebra, Statistics
function main()
# this illustrates the Cholesky decomposition
# create a random p.d. covariance matrix
z = randn(3,3)
V = z'*z # make it pos. def.
println("true covariance matrix")
prettyprint(V)

# take Cholesky decomp.
P = cholesky(V).U
println("P=chol(V): P'P-V")
prettyprint(P'*P - V)

# how to sample from N(0,V)
errs = randn(10000,3)*P
println("sample covariance matrix of 10000 draws from N(0,V)")
prettyprint(cov(errs))

# verify that GLS transformation works
P = inv(cholesky(cov(errs)).U)
errs *= P
println("sample covariance matrix of transformed errors")
prettyprint(cov(errs))
return
end
main()
