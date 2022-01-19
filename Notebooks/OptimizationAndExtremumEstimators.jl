### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# ╔═╡ dbd8ab1c-790d-11ec-05d2-4be402640123
using Econometrics
n = 5
x = [ones(n) 1:n]
β = [10., -1.]
ϵ = randn(n)
y = x*β + ϵ
ols(y,x)

# ╔═╡ 8427ef8a-790e-11ec-256f-25ac76ab14bb
using Statistics
objᵢ = β -> (y-x*β).^2.  	# iᵗʰ obs contrib to sum of squares
obj = β -> mean(objᵢ(β))	# average objective
fminunc(obj, zeros(2))

# ╔═╡ 6da61ee8-790f-11ec-2354-afec042707f3
using Optim
tol = 1e-08
results = Optim.optimize(obj, zeros(2), LBFGS(), 
                            Optim.Options(
                            g_tol = tol,
                            x_tol=tol,
                            f_tol=tol); autodiff=:forward)


# ╔═╡ 30f867ee-7914-11ec-39f2-e1c8568f126b
# the calculations to get I hat
using ForwardDiff
sc =  ForwardDiff.jacobian(objᵢ, βhat)
I = zeros(2,2)
for i = 1:n
	I .+= sc[i,:]*sc[i,:]'
end
I ./= n
# you could also use simply cov(sc), which converges to the same thing.



# ╔═╡ 8e6f67ee-7917-11ec-0f4d-45e17646b9b1
using LinearAlgebra
v∞hat = inv(J)*I*inv(J) # this is the estimate of the limiting var of √N(β-β⁰) 
t = sqrt.(diag(v∞hat/5.))   # to get small sample est. variance, divide by n

# ╔═╡ a7517d2a-790c-11ec-1571-cfd030a3b758
md"# Numeric Optimization"

# ╔═╡ 85b85e6c-790d-11ec-1c31-35cbb8ae8f16
md"## OLS"

# ╔═╡ bcb2d442-790d-11ec-1ca2-931ae3840b28
md"Let's start with a simple problem, where we already know the answers. We will apply numeric optimization methods and extremum estimation theory, and we will be able to verify that they give use the correct results. To start, let's define some data, and do the basic OLS computations"

# ╔═╡ 68b93128-790e-11ec-32e7-d18f5d26cd32
md"Let's see how to get these results using numeric minimization and extremum estimation theory..."

# ╔═╡ 1a9a5c2a-790f-11ec-171c-5b94beb03e82
md"fminunc is just a frontend to the Optim.jl package's LBFGS unconstrained minimizer routine. To call that directly, do"

# ╔═╡ 16cec664-7910-11ec-1a97-7773f65b7042
md"The output has the function value and convergence information. To see the value of the minimizer, do"

# ╔═╡ 47d5f4e4-7910-11ec-3a4b-c5f3378c732c
βhat = results.minimizer

# ╔═╡ f02dacd8-7913-11ec-0e6a-933ca7299a73
md"How do we get the estimated standard errors, using extremum estimation theory? We need I, the covariance of the score contributions, and J, the Hessian matrix."

# ╔═╡ 4be110ee-7917-11ec-03af-01be9183e472
J = ForwardDiff.hessian(obj, βhat)

# ╔═╡ Cell order:
# ╠═a7517d2a-790c-11ec-1571-cfd030a3b758
# ╠═85b85e6c-790d-11ec-1c31-35cbb8ae8f16
# ╠═bcb2d442-790d-11ec-1ca2-931ae3840b28
# ╠═dbd8ab1c-790d-11ec-05d2-4be402640123
# ╠═68b93128-790e-11ec-32e7-d18f5d26cd32
# ╠═8427ef8a-790e-11ec-256f-25ac76ab14bb
# ╠═1a9a5c2a-790f-11ec-171c-5b94beb03e82
# ╠═6da61ee8-790f-11ec-2354-afec042707f3
# ╠═16cec664-7910-11ec-1a97-7773f65b7042
# ╠═47d5f4e4-7910-11ec-3a4b-c5f3378c732c
# ╠═f02dacd8-7913-11ec-0e6a-933ca7299a73
# ╠═30f867ee-7914-11ec-39f2-e1c8568f126b
# ╠═4be110ee-7917-11ec-03af-01be9183e472
# ╠═8e6f67ee-7917-11ec-0f4d-45e17646b9b1
