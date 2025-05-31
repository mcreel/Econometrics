# White and Goldfeld-Quandt tests for heteroscedasticity using
# the Nerlove Cobb-Douglas model
load nerlove.data;

data = data(:,2:6);
data = log(data);
n = rows(data);
y = data(:,1);
x = data(:,2:5);
x = [ones(n,1), x];
k = columns(x);

# create the block diagonal X matrix corresponding to separate coefficients
big_x = zeros(n,5*k);
for i=1:k
	startrow = (i-1)*29+1;
	endrow = i*29;
	startcol =(i-1)*k + 1;
	endcol = i*k;
	big_x(startrow:endrow,startcol:endcol) \
		= big_x(startrow:endrow,startcol:endcol) \
		+ x(startrow:endrow,:);
endfor
x = big_x;

names = char("constant", "output","labor", "fuel", "capital");
names = [names; names; names; names; names]; # copy 5 times

# Now let's try out the model with input prices restricted but constant and output varying
R = eye(5);
R = R(3:5,:);
Z = zeros(3,5);
R = [
	R -R Z Z Z;
	R Z -R Z Z;
	R Z Z -R Z;
	R Z Z Z -R
	];

r = zeros(12,1);

[junk, junk, e] = mc_olsr(y, x, R, r, names);


# WHITE'S TEST
WhiteTest(e, x);

# NOW GOLDFELD-QUANDT (without restrictions) we'll drop the middle 29 obsn
# Also, we put the ESS for subsample of larger firms in the DENOMINATOR
# since the error plot shows that this is the group that has
# the smaller variance
y1 = y(1:58,:);
y3 = y(88:145,:);
x1 = x(1:58,:);
x3 = x(88:145,:);
n = rows(y1);
k = columns(x1);
# the test statistic - simple formula, since the subsamples are equally sized
gq_test = ess(y1,x1) / ess(y3,x3);
gq_pvalue = 1 - fcdf(gq_test,n-k, n-k);
result = [gq_test gq_pvalue];

# print it out
rlabels = char("GQ test");
clabels = char("Value","p-value");
prettyprint(result, rlabels, clabels);

