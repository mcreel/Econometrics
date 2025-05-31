# Estimates 5 separate models for the Nerlove Cobb-Douglas model,
# and does a Chow test for pooled coefficients
using Econometrics, Plots, DelimitedFiles, LinearAlgebra
function main()
data = readdlm("nerlove.data")
data = data[:,2:6]
data = log.(data)
n = size(data,1)
y = data[:,1]
x = data[:,2:end]
x = [ones(n,1) x]
#b, junk, junk, junk = ols(y, x, names=names)
k = size(x,2)

# create the block diagonal X matrix corresponding to separate coefficients
big_x = zeros(n,5*k)
for i=1:k
	startrow = (i-1)*29+1
	endrow = i*29
	startcol =(i-1)*k + 1
	endcol = i*k
	big_x[startrow:endrow,startcol:endcol] += x[startrow:endrow,:]
end
x = big_x

names = ["constant", "output","labor", "fuel", "capital"]
names = [names; names; names; names; names] # copy 5 times

println("Nerlove model: 5 separate regressions for different output levels")
b, junk, junk, junk = ols(y, x, names=names)
output = 1:5
output = output .*5 .+2 .- 5
output = b[output,:]
rts = 1 ./ output

# gset term X11
group = 1:5
group = group
plot(group, rts, label = "", yaxis="RTS", xaxis = "Output Group", show=true)
#savefig("rts.png")

# Chow test
R = Matrix{Float64}(I, 5, 5)
Z = zeros(5,5)
R = [
	R -R Z Z Z;
	R Z -R Z Z;
	R Z Z -R Z;
	R Z Z Z -R]
r = zeros(20,1)

println("Chow test: note that the restricted model")
println("gives the same results as the original model")
ols(y, x, R=R, r=r, names=names)
TestStatistics(y, x, R, r)
return
end
main()

