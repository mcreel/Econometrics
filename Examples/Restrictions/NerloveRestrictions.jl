## Estimates the basic Nerlove Cobb-Douglas model
using Econometrics, DelimitedFiles
cd(@__DIR__)
data = readdlm("nerlove.data")
data = data[:,2:6]
data = log.(data)
n = size(data,1)
y = data[:,1]
x = data[:,2:end]
x = [ones(n,1) x]
names = ["constant", "output", "labor", "fuel", "capital"]
b, junk, junk, junk = ols(y, x, names=names);

## First HOD1
R = [0 0 1 1 1]
r = 1
println()
println("Imposing and testing HOD1")
ols(y, x, R=R, r=r, names=names)
TestStatistics(y, x, R, r);

## Now CRTS
R = [0 1 0 0 0]
r = 1
println()
println("Imposing and testing CRTS")
ols(y, x, R=R, r=r, names=names)
TestStatistics(y, x, R, r);

## Now both, jointly
R = [0 1 0 0 0; 0 0 1 1 1]
r = [1; 1]
println()
println("Imposing and testing HOD1 and CRTS")
ols(y, x, R=R, r=r, names=names)
TestStatistics(y, x, R, r);
