#
# 	The dep. vbls, in corresponding column
# 	Office based doctor visits	    1
# 	Outpatient doctor visits		2
# 	Emergency room visits			3
# 	Inpatient visits				4
# 	Dental visits					5
# 	Prescriptions					6
#

1;

which_dep = 1;
if (which_dep == 1) printf("\nOBDV\n"); endif
if (which_dep == 2) printf("\nOPV");endif
if (which_dep == 3) printf("\nIPV");endif
if (which_dep == 4) printf("\nERV");endif
if (which_dep == 5) printf("\nDV");endif
if (which_dep == 6) printf("\nPRESCR");endif

load meps1996.data;
y = data(:,which_dep);
x = data(:,7:12);
n = rows(x);
x = [ones(n,1) x];
[x scale] = scale_data(x); # scale the data

names = char("constant","pub. ins.","priv. ins.", "sex", "age","edu","inc");

function m = qml(theta, data, momentargs)
	kx = momentargs{1};
	data = data(:,1:kx+1);
	[junk, m] = Poisson(theta, data); # moments are Poisson scores
endfunction

function m = iv(theta, data, momentargs)
	kx = momentargs{1};
	y = data(:,1);
	x = data(:,2:kx+1);
	w = data(:,kx+2:columns(data));
	e = y .* exp(-x*theta) - 1; % preferred, it's less heteroscedastic
	%e = y - exp(x*theta);     % not preferred, but legitimate
	m = diag(e)*w;
endfunction



# define instruments
k = columns(x);
w = [x(:,2) , x(:,4:k)];
z = x(:,4:k);
zz = diag(x(:,2))*z;
w = [w, zz];
[w, junk, junk2] = st_norm(w);
w = [ones(n,1), w];
kx = columns(x);
kw = columns(w);
momentargs = {kx, kw};
data = [y x w];
control = {-1};

# IV estimator
kk = columns(w);
theta = zeros(k,1);
moments = "iv";
weight = eye(kk);
for i = 1:10
	theta = gmm_estimate(theta, data, weight, moments, momentargs);
	# optimal weight, to get correct varcov
	m = feval(moments, theta, data, momentargs);
	weight = inv(cov(m));
endfor
title = "IV";
gmm_results(theta, data, weight, moments, momentargs, names, title, scale);

# Poisson QML estimator, estimated in GMM form
k = columns(x);
moments = "qml";
weight = eye(k);
theta = gmm_estimate(theta, data, weight, moments, momentargs);
# optimal weight, to get correct varcov
m = feval(moments, theta, data, momentargs);
weight = inv(cov(m));
title = "Poisson QML";
gmm_results(theta, data, weight, moments, momentargs, names, title, scale);


