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

# Plot the fitted line
#plot(x,PopRegLine, label = "population regression line")
xlabel= "X"
scatter(x, y, label = "data")
#savefig("TypicalData.png")
#gui()
