using Plots
# plots the contours of the OLS objective function,
# with the true and estimated coefficients
# generate data
n = 20
x = [ones(n) randn(n)]
beta = rand(2)
y = x*beta+randn(n)
# ols estimate and objective function
betahat = x\y
s = (a,b)-> (1/n)*sum((y - x*[a,b]).^2)
# plot the contours and the points
closeall()
b1 = range(-1.0,stop=1.0,length=100) 
b2 = range(-1.0,stop=1.0,length=100)
p = contour(b1,b2,(b1,b2)->s(b1,b2),fill=true, c=:viridis)
scatter!([beta[1]],[beta[2]], markersize=10, label="true")
scatter!([betahat[1]],[betahat[2]], markersize=10, label="estimated")
