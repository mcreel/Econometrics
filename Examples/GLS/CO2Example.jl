using Econometrics, DelimitedFiles, Plots
cd(@__DIR__)
CO2 = readdlm("../Data/CO2.data")
plot(CO2)
n = size(CO2,1)
trend = 1:n
x = [ones(n) trend]
b, junk, e, junk, junk= ols(CO2,x)
fit = x*b
p = plot(trend[end-36:end], e[end-36:end], legend=false)
title!("OLS residuals, last 3 years of data")
display(p)
#savefig("CO2Residuals.png")
