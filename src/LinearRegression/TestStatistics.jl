# TESTS - the classical test statistics for linear restrictions
# 
# written by M. Creel, March 10, 1998.  michael.creel@uab.es
# converted to Octave 13 Jan. 2003
# converted to Julia 25 Aug. 2017
# 
# purpose: calculates qF, Wald, Score and Likelihood Ratio tests for
# linear model						 y=XB+e
# 									 e~N(0,sig^2*I_n)
# subject to linear restrictions  	 RB=r
# 
# format: [qF, W, LR, S] = TestStatistics(y,x,R,r)
# 
# inputs:
# 		 y: nx1 dependent variable
# 		 x: nxk regressor matrix
# 		 R: R above, a qxk matrix
# 		 r: r above, a qx1 vector
# 
# returns: 
#        qF: the qF statistic
# 		 W: the Wald statistic
# 		 S: the score statistic
# 		 LR: the likelihood ratio statistic
using StatsFuns
function TestStatistics(y, x, R, r; silent=false)
    n,k = size(x)
    q = size(R,1)
	b = x\y
	xx_inv = inv(x'*x)
	P_inv = inv(R*xx_inv*R')
	b_r = b .- xx_inv*R'*P_inv*(R*b.-r)
    e = y - x*b
    ess = (e'*e)[1]
    e_r = y - x*b_r
    ess_r = (e_r' * e_r)[1]
    sigsqhat_ols = ess/(n-k)
    sigsqhat_mle = ess/(n)
    sigsqhat_mle_r = ess_r/(n)
    # F-test
    F = (ess_r-ess)/(q*sigsqhat_ols)
    # Wald test (uses unrestricted model's est. of sig^2
    W = (R*b.-r)'*P_inv*(R*b.-r)/sigsqhat_mle
    # Score test (uses restricted model's est. of sig^2 
    P_x = x * xx_inv * x'
    S = e_r' * P_x * e_r/(sigsqhat_mle_r)
    # LR test
    lnl = -n/2*log(2*pi) - n/2*log(sigsqhat_mle) - ess/(2.0*sigsqhat_mle)
    lnl_r = -n/2*log(2*pi) - n/2*log(sigsqhat_mle_r) - ess_r/(2.0*sigsqhat_mle_r)
    LR = 2.0*(lnl-lnl_r)
	if !silent
        tests = [q*F[1], W[1], LR[1], S[1]]
	    pvalues = chisqccdf.(tests,q)
        tests = [tests pvalues]
	    TESTS = ["qF","Wald","LR","Score"]
	    labels = ["Value","p-value"]
        prettyprint(tests, labels, TESTS)
    end
    return F, W, LR, S
end
