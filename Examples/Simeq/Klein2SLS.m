# Estimates the Klein consumption equation by
# 2SLS, plots residuals, and does
# Breusch-Godfrey test

load klein.data;
data = klein;


# construct missing lags, and drop first row that has missing data
profits = data(:,3);
output = data(:,7);
data = [data lag(profits,1) lag(output,1)];
data = data(2:rows(data),:);

n = rows(data);

# define instruments
exogs = [1, 6, 8, 9, 10, 11, 12];
exogs = data(:,exogs);
exogs = [ones(n,1) exogs];

# CONSUMPTION
printf("CONSUMPTION EQUATION\n");
# define variables in consumption equation
y = data(:,2);
profits = data(:,3);
lagprofits = data(:,11);
wp = data(:,4);
wg = data(:,8);
wages = wp + wg;
# regressors in consumption equation
x = [profits lagprofits wages];
x = [ones(n,1) x];
# get the 2SLS instruments
b = ols(x, exogs);
xhat = exogs*b;
# 2SLS estimation
names = char("Constant", "Profits", "Lagged Profits", "Wages");
[junk, junk, e] = mc_2sls(y, x, xhat, names);


# INVESTMENT
printf("INVESTMENT EQUATION\n");
# define variables in investment equation
y = data(:,5);
profits = data(:,3);
lagprofits = data(:,11);
lagcapital = data(:,6);
# regressors in  investment equation
x = [profits lagprofits lagcapital];
x = [ones(n,1) x];
# get the 2SLS instruments
b = ols(x, exogs);
xhat = exogs*b;
# 2SLS estimation
names = char("Constant", "Profits", "Lagged Profits", "Lagged Capital");
[junk, junk, e] = mc_2sls(y, x, xhat, names);


# WAGES
printf("WAGES EQUATION\n");
# define variables in wages equation
y = data(:,4);
year = data(:,1)-1931;
output = data(:,7);
lagoutput = data(:,12);
# regressors in wages equation
x = [output lagoutput year];
x = [ones(n,1) x];
# get the 2SLS instruments
b = ols(x, exogs);
xhat = exogs*b;
# 2SLS estimation
names = char("Constant", "Output", "Lagged Output", "Trend");
[junk, junk, e] = mc_2sls(y, x, xhat, names);
