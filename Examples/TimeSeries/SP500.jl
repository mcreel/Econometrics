using CSV, DataFrames, Econometrics, StatsModels, Statistics
sp500 = CSV.read("../Data/sp500.csv", DataFrame)
y = (sp500.rets).^2
n = size(y,1)
x = [ones(n) lags(y,2)]
y = y[3:end]
x = x[3:end,:]
b, varb, e, junk, junk = ols(y,x,vc="nw")
e = [e lags(e,8)]
cor(e[9:end,:])