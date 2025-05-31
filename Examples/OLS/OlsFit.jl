# this generates data according to classical model
# and shows the OLS fit
using Plots
n = 20
x = 1:n

X = [ones(n,1) x]
beta = [10, -1]
e = 6.0*randn(n)
PopRegLine = X*beta
y = PopRegLine + e

# the OLS coefficients
b = inv(X'*X)*X'*y
fit = X*b

# Plot the fitted line
#plot(x,PopRegLine, label = "population regression line")
plot(x, fit, label = "OLS fitted line")
xlabel= "X"
scatter!(x, y, label = "data")
#savefig("OlsFit.png")
#gui()
