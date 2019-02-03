# this performs a kernel regression estimation of E(y^2(t)) conditional on y(t-1) and y(t-2)
# using yen/dollar and euro/dollar spot rates and then plots the surfaces


points = 30;
ww = "";
close all;


# euro/dollar
load rates;
data = rates(:,1);
maxlag = 2;
data = data .^2;				# square it to look at volatility
data = [data lags(data,maxlag)];		# make lags
data = data(maxlag+1:end,:);			# drop missings
data = st_norm(data);				# mean is now zero
y = data(:,1);					# dependent variable
x = data(:,2:maxlag+1); 			# the two lags are conditioning vbls
x = 6*normcdf(x) - 3;				# compactify - zero maps to zero
# evaluation points are on a grid for plotting
a = -3:0.2:3;
a = a';
gridsize = rows(a);
[grid_x, grid_y] = meshgrid(a, a);
eval_points = [vec(grid_x) vec(grid_y)];
# the kernel fit and the plot
fit = kernel_regression(eval_points, y, x, ww);
fit = reshape(fit, gridsize, gridsize);
figure();
surf(grid_x, grid_y, fit);
title("Euro/Dollar spot rate: y^2(t)");
xlabel("y^2(t-1)");
ylabel("y^2(t-2)");
%print("eurodollar.png", "-dpng"); # let's not overwrite that file!


# yen/dollar
load rates;
data = rates(:,2);
maxlag = 2;
data = data .^2;				# square it to look at volatility
data = [data lags(data,maxlag)];		# make lags
data = data(maxlag+1:rows(data),:);		# drop missings
data = st_norm(data);				# mean is now zero
y = data(:,1);					# dependent variable
x = data(:,2:maxlag+1); 			# the two lags are conditioning vbls
x = 6*normcdf(x) - 3;				# compactify - zero maps to zero
# evaluation points are on a grid for plotting
a = -3:0.2:3;
a = a';
gridsize = rows(a);
[grid_x, grid_y] = meshgrid(a, a);
eval_points = [vec(grid_x) vec(grid_y)];
# the kernel fit and the plot
fit = kernel_regression(eval_points, y, x, ww);
fit = reshape(fit, gridsize, gridsize);
figure();
surf(grid_x, grid_y, fit);
title("Yen/Dollar spot rate: y^2(t)");
xlabel("y^2(t-1)");
ylabel("y^2(t-2)");
%print("yendollar.png", "-dpng"); # let's not overwrite that file!
