# this program illustrates sample selection bias associated with
# dropping observations for which the dep. vbl. is <0. The resulting figure is
# saved in sampsel.ps
using Plots
pyplot()
n = 100
sig = 5
x = 10*rand(n,1)
y = x+sig*randn(n,1)
# this drops the rows for which y < 0
z = y.>0.0
yin = y[z]  
xin = x[z]
# regression using selected sample
xin = [ones(size(xin)) xin]
b = xin\yin
# plots
xx=(0:0.05:10)
yy=xx
yhat = b[1]+b[2]*xx
scatter(x,y, label = "data")
plot!(xx,yy,label = "population regression line")
plot!(xx,yhat,label = "fitted regression line")
gui()
savefig("sampsel.png")


