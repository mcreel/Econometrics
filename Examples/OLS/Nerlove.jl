# Estimates the basic Nerlove Cobb-Douglas model
using DelimitedFiles
function main()
data = readdlm("nerlove.data")
data = data[:,2:6]
data = log.(data)
n = size(data,1)
y = data[:,1]
x = data[:,2:end]
x = [ones(n,1) x]
names = ["constant", "output", "labor", "fuel", "capital"]
b, junk, junk, junk = ols(y, x, names=names)
return
end
main()
