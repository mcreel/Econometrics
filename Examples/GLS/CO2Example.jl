using Econometrics, DelimitedFiles, Plots
CO2 = readdlm("../Data/CO2.data")
plot(CO2)
n = size(CO2,1)
trend = 1:n
x = [ones(n) trend]
b, junk, e, junk, junk= ols(CO2,x)
fit = x*b
plot(trend[end-36:end], e[end-36:end], legend=false)
title!("OLS residuals, last 3 years of data")
savefig("CO2Residuals.png")
gui()
