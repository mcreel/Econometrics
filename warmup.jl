# this script will precompile some commonly used functions
using Econometrics, StatsPlots, CSV, DataFrames
ols()
mleresults()
gmmresults()
samin()
fminunc()
fmincon()
mcmc()
npreg()
dstats(rand(100,3))
plot(1:10)
histogram(rand(100))
scatter(1:10)
d = DataFrame(rand(2,2), :auto)
CSV.write(tempdir()*"/junk.csv", d)
d = CSV.read(tempdir()*"/junk.csv", DataFrame)
clc()

