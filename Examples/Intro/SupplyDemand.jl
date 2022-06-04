## this generates data from the simple supply-demand model
#  q = α1 + α2*p + α3*m + e1
#  q = β1 + β2*p +        e2
## setup
using Econometrics, Plots, DelimitedFiles
# number of obsn
n = 500
# model parameters
α1 = 100.0
α2 = -1.0
α3 = 1.0
β1 = 20.0
β2 = 1.0

## make variables and shocks
# exog var income
m = randn(n,5)
m = sum(m.*m,dims=2) # chi square df=5
# structural shocks
ϵ1 = 1.0*randn(n)
ϵ2 = 1.0*randn(n)
# rf shocks
ν1 = (β2*ϵ1-α2*ϵ2)/(β2-α2)
ν2 = (ϵ1-ϵ2)/(β2-α2)
# rf coefs
π11 = (β2*α1-α2*β1)/(β2-α2)
π21 = β2*α3/(β2-α2)
π12 = (α1-β1)/(β2-α2)
π22 = α3/(β2-α2)

##  generate endog variables computed from rf
q = π11 .+ π21*m + ν1
p = π12 .+ π22*m + ν2
## write and plot data
data = [q p m]
data = sortbyc(data,3)
nn = Int(n/4)
p = scatter(data[1:nn,2], data[1:nn,1], markersize=8, xlabel="price", ylabel="quantity", legend=false)
for i = 2:4
    scatter!(data[(i-1)*nn+1:i*nn,2], data[(i-1)*nn:i*nn,1], markersize=8, show=true)
end
display(p)
#writedlm("data.txt", data)
