#!/usr/bin/env python
from scipy import *
import scipy.io.array_import
import scipy.stats.kde
from pylab import *

filename=('rates')
data = scipy.io.array_import.read_array(filename)


# dollar-euro
# the histogram of the data
fig=figure()
x = data[:,0]
mu = mean(x)
sigma = std(x)
n, bins, patches = hist(x, 50, normed=1)
setp(patches, 'facecolor', 'g', 'alpha', 0.75)
# add a 'best fit' line
y = normpdf( bins, mu, sigma)
l = plot(bins, y, 'r--')
setp(l, 'linewidth', 1)
xlabel('rate')
ylabel('Density')
title(r'$\rm{Growth\ rate\ of\ dollar/euro\ exchange\ rate}\ $')
grid(True)

# time series plot
fig=figure()
xlabel('time')
ylabel('rate')
title(r'$\rm{Growth\ rate\ of\ dollar/euro\ exchange\ rate}\ $')
plot(x)



# dollar-yen
# the histogram of the data
fig=figure()
x = data[:,1]
mu = mean(x)
sigma = std(x)
n, bins, patches = hist(x, 50, normed=1)
setp(patches, 'facecolor', 'g', 'alpha', 0.75)

# add a 'best fit' line
y = normpdf( bins, mu, sigma)
l = plot(bins, y, 'r--')
setp(l, 'linewidth', 1)

xlabel('rate')
ylabel('Density')
title(r'$\rm{Growth\ rate\ of\ dollar/yen\ exchange\ rate}\ $')
#axis([40, 160, 0, 0.03])
grid(True)


# time series plot
fig=figure()
xlabel('time')
ylabel('rate')
title(r'$\rm{Growth\ rate\ of\ dollar/yen\ exchange\ rate}\ $')
plot(x)


show()

