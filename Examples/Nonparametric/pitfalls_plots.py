#!/usr/bin/env python
from pylab import *


## linear fit
fig=figure()
x = arange(0.0, 2*pi, 0.01)
true = 1.0 + 3.0*x/(2.0*pi) - (x/(2.0*pi))**2.0
approx = 7.0/6.0 + x/pi
plot(x, approx, 'r', label='approx')
plot(x, true, 'b', label='true')
xx = concatenate( (x,x[::-1]) )
y = concatenate( (true,approx[::-1]) )
p = fill(xx, y, facecolor='g', alpha=0.5)
legend(loc='best')
xlabel('x')
#fig.savefig('linearfit.eps')


## linear elasticity
fig=figure()
true_derivative = 3.0/(2.0*pi) - 2.0*x/(2.0*pi)/(2.0*pi)
true_elasticity = x*true_derivative/true
approx_derivative = 1.0/pi
approx_elasticity = x*approx_derivative/approx
plot(x, approx_elasticity, 'r', label='approx')
plot(x, true_elasticity, 'b', label='true')
xx = concatenate( (x,x[::-1]) )
y = concatenate( (true_elasticity,approx_elasticity[::-1]) )
p = fill(xx, y, facecolor='g', alpha=0.5)
legend(loc='best')
xlabel('x')
#fig.savefig('linearelasticity.eps')


## fourier fit
fig=figure()
true = 1.0 + 3.0*x/(2.0*pi) - (x/(2.0*pi))**2.0
approx = 7.0/6.0 + x/pi -cos(x)/(pi**2) - cos(2*x)/(4*pi**2)
plot(x, approx, 'r', label='approx')
plot(x, true, 'b', label='true')
xx = concatenate( (x,x[::-1]) )
y = concatenate( (true,approx[::-1]) )
p = fill(xx, y, facecolor='g', alpha=0.5)
legend(loc='best')
xlabel('x')
#fig.savefig('fourierfit.eps')


## fourier elasticity
fig=figure()
true_derivative = 3.0/(2.0*pi) - 2.0*x/(2.0*pi)/(2.0*pi)
true_elasticity = x*true_derivative/true
approx_derivative = 1.0/pi + sin(x)/(pi**2) + 2*sin(2*x)/(4*pi**2)
approx_elasticity = x*approx_derivative/approx
plot(x, approx_elasticity, 'r', label='approx')
plot(x, true_elasticity, 'b', label='true')
xx = concatenate( (x,x[::-1]) )
y = concatenate( (true_elasticity,approx_elasticity[::-1]) )
p = fill(xx, y, facecolor='g', alpha=0.5)
legend(loc='best')
xlabel('x')
#fig.savefig('fourierelasticity.eps')


show()
