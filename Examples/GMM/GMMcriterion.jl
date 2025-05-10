# this shows how the GMM criterions combines information from two
# moment conditions. 
using Plots, LinearAlgebra

function main(weight_on_first)
# the two separate moment conditions
m1(θ) = 2 + 3θ
m2(θ) = 5 -3θ + 0.1*θ^2 

# the GMM criterions
m(θ) = [m1(θ), m2(θ)]
W = 0.05*[weight_on_first 0. ; 0. 1. - weight_on_first]
s(θ) = dot(m(θ),W,m(θ))[1]

# plot it
plot([m1,m2,s], labels=["m¹ₙ(θ)" "m²ₙ(θ)" "sₙ(θ) "], title="Two moment conditions, and the GMM criterion")

# savefig("gmmcriterion.png")
end

main(0.0)  # set the weight on first component between zero and one