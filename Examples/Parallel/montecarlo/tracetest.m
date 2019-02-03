# Simulates trace test for cointegration 
# Can be used to tabulate empirical distribution
#
# Ref. Doornik et. al. (2002) 
# "Computationally-intensive Econometrics using
# a Distributed Matrix-programming Language",
# Philosophical Transactions of the Royal Society of London,
# Series A, 360, 1245-1266. 
function test_stat = tracetest(args)
	t = args{1};
	n = args{2};
	e = randn(t,n);
	p = inv(chol(e'*e/t));
	e = e*p;
	s = lag(cumsum(e),1);
	s(1,:) = s(1,:) - s(1,:); # lag fills with 1, test needs 0
	fac = e'*s;
	ev = eig(fac*inv(s'*s)*fac'/t);
	test_stat = -t*sum(log(1 - ev/t));
endfunction 
