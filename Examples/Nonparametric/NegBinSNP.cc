// Copyright (C) 2004   Michael Creel   <michael.creel@uab.es>
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// This is the NegBinSNP loglikelihood function

#include <octave/oct.h>
#include <octave/Cell.h>
#include <octave/parse.h>
#include <float.h>


// define argument checks
static bool
any_bad_argument(const octave_value_list& args)
{
	if (!args(0).is_matrix_type())
    	{
        	error(": NegBinSNP: first argument vector theta (parameters)");
        	return true;
    	}
	if (!args(1).is_matrix_type())
	{
        	error("NegBinSNP: second argument must be a matrix of data");
        	return true;
	}
	Cell otherargs (args(2).cell_value());
	if (!(otherargs.length() == 1))
	{
	 	error("NegBinSNP: the third argument must be cell");
 		return true;
	}
	if (!otherargs(0).is_real_scalar())
	{
        	error("NegBinSNP: third arg must cell containing negbin_type (1/2)");
        	return true;
	}
    	return false;
}

DEFUN_DLD(NegBinSNP, args, ,"NebBinSNP likelihood function")
{

	 // check the arguments
 	if (any_bad_argument(args)) return octave_value_list();

	ColumnVector theta = args(0).column_vector_value();
	Matrix data (args(1).matrix_value());

	const int n = data.rows();
	const int k = data.columns() - 1;
	ColumnVector y = data.column(0);

	Matrix x = data.extract_n(0, 1, n, k);

	Cell otherargs (args(2).cell_value());
	int nb_type (otherargs(0).int_value());
	octave_value_list f_return;

	int i, j;
	ColumnVector m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11;
	double t1, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16;
	double t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27;
	double t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39;
	double t40, t41, t42, t43, t44;

	const int kk = theta.rows();
	const int order = kk - k - 1; // that "1" is for alpha

	if (order > 4)
	{
        error("NegBinSNP: polynomial order is limited to 4 at most, at present");
        return octave_value_list();
	}


	double alpha = DBL_EPSILON + exp(theta(k));
	ColumnVector beta(k);
	ColumnVector gam(order + 1); // 2nd degree has 3 coefficients, e.g.
	ColumnVector lam(n);
	ColumnVector psi(n);
	ColumnVector lampsi(n);
	ColumnVector logdensity(n);
	ColumnVector xbeta;
	ColumnVector normfactor(n);
	ColumnVector prediction(n);
	ColumnVector poly(n);

	// how many moments are needed?
	// two times order for density,
	// two times order  plus 1 for prediction
	Matrix moments(n,2*order+2);
	moments.fill(1.0);

	// extract beta as leading k elements of theta
	for(i = 0; i < k; i++)
	{
		beta(i) = theta(i);
	}

	gam(0) = 1; // normalization
	// extract gam as last kk - k  - 1 elements of theta
	for(i = 0; i < order; i++)
	{
		gam(i+1) = theta(k+1+i) / pow(10.0, ((double) i) +1.0);
	}

	xbeta = x*beta;

	for (i=0; i < n; i++)
	{
		lam(i) = DBL_EPSILON + exp(xbeta(i));
	}

	if (nb_type==1)
	{
		psi = lam / (alpha);
	}
	if (nb_type==2)
	{
		psi = psi.fill(DBL_EPSILON + 1.0/alpha);
	}


	// top-bound psi to bulletproof
	for (i=0; i < n; i++)
	{

		psi(i) = (psi(i) < 10000) * psi(i) + (psi(i) > 10000)*10000;
	}


	// moments
	if (order==1)
	{
		for(i = 0; i < n; i++)
		{
			moments(i,1) = lam(i);

    			t11 = lam(i) + psi(i) ;
    			t12 = psi(i)/t11 ;
    			t13 = pow(t12, psi(i)) ;
    			t14 = -psi(i) - 1.0 ;
    			t15 = -(lam(i)/t11) + 1.0 ;
    			t16 = log(t15) ;
    			moments(i,2) = (t13*lam(i)*psi(i)*exp(t14*t16))/t11 - pow(t11, -2.0)*t13*t14*(lam(i)*lam(i)) \
					*psi(i)*exp(t16*(-psi(i) - 2.0)) ;


    			t15 = lam(i) + psi(i) ;
    			t16 = psi(i)/t15 ;
    			t17 = pow(t16, psi(i)) ;
    			t18 = -psi(i) - 1.0 ;
    			t24 = lam(i)/t15 ;
    			t19 = -t24 + 1.0 ;
    			t20 = log(t19) ;
    			t21 = lam(i)*lam(i) ;
    			t22 = t15*t15 ;
    			t23 = -psi(i) - 2.0 ;
    			moments(i,3) = (t17*lam(i)*psi(i)*exp(t20*t18))/t15 - 3.0*t21*pow(t15, -2.0)*t17*t18 \
				* psi(i)*exp(t20*t23) + (t21*t23*t17*t18*lam(i)*psi(i)*exp(t20*(-psi(i) - 3.0)))/(t22*t15) ;
		}
	}

	if (order==2)
	{
		for(i = 0; i < n; i++)
		{
			moments(i,1) = lam(i);

    			t11 = lam(i) + psi(i) ;
    			t12 = psi(i)/t11 ;
    			t13 = pow(t12, psi(i)) ;
    			t14 = -psi(i) - 1.0 ;
    			t15 = -(lam(i)/t11) + 1.0 ;
    			t16 = log(t15) ;
    			moments(i,2) = (t13*lam(i)*psi(i)*exp(t14*t16))/t11 - pow(t11, -2.0)*t13*t14*(lam(i)*lam(i)) \
					*psi(i)*exp(t16*(-psi(i) - 2.0)) ;


    			t15 = lam(i) + psi(i) ;
    			t16 = psi(i)/t15 ;
    			t17 = pow(t16, psi(i)) ;
    			t18 = -psi(i) - 1.0 ;
    			t24 = lam(i)/t15 ;
    			t19 = -t24 + 1.0 ;
    			t20 = log(t19) ;
    			t21 = lam(i)*lam(i) ;
    			t22 = t15*t15 ;
    			t23 = -psi(i) - 2.0 ;
    			moments(i,3) = (t17*lam(i)*psi(i)*exp(t20*t18))/t15 - 3.0*t21*pow(t15, -2.0)*t17*t18 \
				* psi(i)*exp(t20*t23) + (t21*t23*t17*t18*lam(i)*psi(i)*exp(t20*(-psi(i) - 3.0)))/(t22*t15) ;

 	  		t16 = lam(i) + psi(i) ;
    			t17 = psi(i)/t16 ;
    			t18 = pow(t17, psi(i)) ;
    			t19 = -psi(i) - 1.0 ;
    			t25 = lam(i)/t16 ;
    			t20 = -t25 + 1.0 ;
    			t21 = log(t20) ;
    			t22 = lam(i)*lam(i) ;
    			t23 = t16*t16 ;
    			t24 = -psi(i) - 2.0 ;
    			t26 = -psi(i) - 3.0 ;
    			moments(i,4) = (t18*lam(i)*psi(i)*exp(t21*t19))/t16 - 7.0*t22*pow(t16, -2.0)*t18*t19 \
				*psi(i)*exp(t21*t24) + (6.0*t22*t24*t18*t19*lam(i)*psi(i)*exp(t21*t26))/(t23*t16) \
				-(t22*t22)*pow(t23, -2.0)*t24*t26*t18*t19*psi(i)*exp(t21*(-psi(i) - 4.0)) ;

    			t19 = lam(i) + psi(i) ;
    			t20 = psi(i)/t19 ;
    			t21 = pow(t20, psi(i)) ;
    			t22 = -psi(i) - 1.0 ;
    			t28 = lam(i)/t19 ;
    			t23 = -t28 + 1.0 ;
    			t24 = log(t23) ;
    			t25 = lam(i)*lam(i) ;
    			t26 = t19*t19 ;
    			t27 = -psi(i) - 2.0 ;
    			t29 = -psi(i) - 3.0 ;
    			t30 = t25*t25 ;
    			t31 = t26*t26 ;
    			t32 = -psi(i) - 4.0 ;
    			moments(i,5) = (t21*lam(i)*psi(i)*exp(t22*t24))/t19 - 15.0*t21*t22*t25*pow(t19, -2.0)\
					*psi(i)*exp(t24*t27) - 10.0*t21*t30*t22*pow(t26, -2.0)*t27*t29*psi(i) \
					*exp(t32*t24) + (25.0*t21*t22*t25*t27*lam(i)*psi(i)*exp(t24*t29))/(t26*t19) \
					+ (t21*t30*t22*t32*t27*t29*lam(i)*psi(i)*exp(t24*(-psi(i) - 5.0)))/(t31*t19) ;
		}
	}

	if (order==3)
	{
		for(i = 0; i < n; i++)
		{
			moments(i,1) = lam(i);

    			t11 = lam(i) + psi(i) ;
    			t12 = psi(i)/t11 ;
    			t13 = pow(t12, psi(i)) ;
    			t14 = -psi(i) - 1.0 ;
    			t15 = -(lam(i)/t11) + 1.0 ;
    			t16 = log(t15) ;
    			moments(i,2) = (t13*lam(i)*psi(i)*exp(t14*t16))/t11 - pow(t11, -2.0)*t13*t14*(lam(i)*lam(i)) \
					*psi(i)*exp(t16*(-psi(i) - 2.0)) ;


    			t15 = lam(i) + psi(i) ;
    			t16 = psi(i)/t15 ;
    			t17 = pow(t16, psi(i)) ;
    			t18 = -psi(i) - 1.0 ;
    			t24 = lam(i)/t15 ;
    			t19 = -t24 + 1.0 ;
    			t20 = log(t19) ;
    			t21 = lam(i)*lam(i) ;
    			t22 = t15*t15 ;
    			t23 = -psi(i) - 2.0 ;
    			moments(i,3) = (t17*lam(i)*psi(i)*exp(t20*t18))/t15 - 3.0*t21*pow(t15, -2.0)*t17*t18 \
				* psi(i)*exp(t20*t23) + (t21*t23*t17*t18*lam(i)*psi(i)*exp(t20*(-psi(i) - 3.0)))/(t22*t15) ;

 	  			t16 = lam(i) + psi(i) ;
    			t17 = psi(i)/t16 ;
    			t18 = pow(t17, psi(i)) ;
    			t19 = -psi(i) - 1.0 ;
    			t25 = lam(i)/t16 ;
    			t20 = -t25 + 1.0 ;
    			t21 = log(t20) ;
    			t22 = lam(i)*lam(i) ;
    			t23 = t16*t16 ;
    			t24 = -psi(i) - 2.0 ;
    			t26 = -psi(i) - 3.0 ;
    			moments(i,4) = (t18*lam(i)*psi(i)*exp(t21*t19))/t16 - 7.0*t22*pow(t16, -2.0)*t18*t19 \
				*psi(i)*exp(t21*t24) + (6.0*t22*t24*t18*t19*lam(i)*psi(i)*exp(t21*t26))/(t23*t16) \
				-(t22*t22)*pow(t23, -2.0)*t24*t26*t18*t19*psi(i)*exp(t21*(-psi(i) - 4.0)) ;

    			t19 = lam(i) + psi(i) ;
    			t20 = psi(i)/t19 ;
    			t21 = pow(t20, psi(i)) ;
    			t22 = -psi(i) - 1.0 ;
    			t28 = lam(i)/t19 ;
    			t23 = -t28 + 1.0 ;
    			t24 = log(t23) ;
    			t25 = lam(i)*lam(i) ;
    			t26 = t19*t19 ;
    			t27 = -psi(i) - 2.0 ;
    			t29 = -psi(i) - 3.0 ;
    			t30 = t25*t25 ;
    			t31 = t26*t26 ;
    			t32 = -psi(i) - 4.0 ;
    			moments(i,5) = (t21*lam(i)*psi(i)*exp(t22*t24))/t19 - 15.0*t21*t22*t25*pow(t19, -2.0)\
					*psi(i)*exp(t24*t27) - 10.0*t21*t30*t22*pow(t26, -2.0)*t27*t29*psi(i) \
					*exp(t32*t24) + (25.0*t21*t22*t25*t27*lam(i)*psi(i)*exp(t24*t29))/(t26*t19) \
					+ (t21*t30*t22*t32*t27*t29*lam(i)*psi(i)*exp(t24*(-psi(i) - 5.0)))/(t31*t19) ;

    			t20 = lam(i) + psi(i) ;
    			t21 = psi(i)/t20 ;
    			t22 = pow(t21, psi(i)) ;
    			t23 = -psi(i) - 1.0 ;
    			t29 = lam(i)/t20 ;
    			t24 = -t29 + 1.0 ;
    			t25 = log(t24) ;
    			t26 = lam(i)*lam(i) ;
    			t27 = t20*t20 ;
    			t28 = -psi(i) - 2.0 ;
    			t30 = -psi(i) - 3.0 ;
    			t31 = t26*t26 ;
    			t32 = t27*t27 ;
    			t33 = -psi(i) - 4.0 ;
    			t34 = -psi(i) - 5.0 ;
    			moments(i,6) = (t22*lam(i)*psi(i)*exp(t23*t25))/t20 - 31.0*pow(t20, -2.0)*t22*t23*t26\
					*psi(i)*exp(t25*t28) - 65.0*t30*t22*t31*t23*pow(t27, -2.0)*t28*psi(i)*exp(t33 \
					*t25) + (90.0*t22*t23*t26*t28*lam(i)*psi(i)*exp(t30*t25))/(t20*t27) \
					+ (15.0*t30*t22*t31*t23*t33*t28*lam(i)*psi(i)*exp(t25*t34))/(t20*t32) \
					- (t30*t22*t31*t23*t33*t34*t26*t28*psi(i)*exp(t25*(-psi(i) - 6.0)))/(t32*t27) ;

    			t21 = lam(i) + psi(i) ;
    			t22 = psi(i)/t21 ;
    			t23 = pow(t22, psi(i)) ;
    			t24 = -psi(i) - 1.0 ;
    			t30 = lam(i)/t21 ;
    			t25 = -t30 + 1.0 ;
    			t26 = log(t25) ;
    			t27 = lam(i)*lam(i) ;
    			t28 = t21*t21 ;
    			t29 = -psi(i) - 2.0 ;
    			t31 = -psi(i) - 3.0 ;
    			t32 = t27*t27 ;
    			t33 = t28*t28 ;
    			t34 = -psi(i) - 4.0 ;
    			t35 = -psi(i) - 5.0 ;
    			t36 = -psi(i) - 6.0 ;
    			moments(i,7) = (t23*lam(i)*psi(i)*exp(t24*t26))/t21 - 63.0*pow(t21, -2.0)*t23*t24*t27 \
					*psi(i)*exp(t26*t29) - 350.0*t31*t23*t32*t24*pow(t28, -2.0)*t29*psi(i) \
					*exp(t34*t26) + (301.0*t23*t24*t27*t29*lam(i)*psi(i)*exp(t31*t26))/(t21*t28) \
					+ (140.0*t31*t23*t32*t24*t34*t29*lam(i)*psi(i)*exp(t26*t35))/(t21*t33) \
					- (21.0*t31*t23*t32*t24*t34*t35*t27*t29*psi(i)*exp(t26*t36))/(t33*t28) \
					+ (t31*t23*t32*t24*t34*t35*t27*t36*t29*lam(i)*psi(i)*exp(t26*(-psi(i) - 7.0)))/(t21*t33*t28) ;
		}
	}


	if (order==4)
	{
		for(i = 0; i < n; i++)
		{
			moments(i,1) = lam(i);

    			t11 = lam(i) + psi(i) ;
    			t12 = psi(i)/t11 ;
    			t13 = pow(t12, psi(i)) ;
    			t14 = -psi(i) - 1.0 ;
    			t15 = -(lam(i)/t11) + 1.0 ;
    			t16 = log(t15) ;
    			moments(i,2) = (t13*lam(i)*psi(i)*exp(t14*t16))/t11 - pow(t11, -2.0)*t13*t14*(lam(i)*lam(i)) \
					*psi(i)*exp(t16*(-psi(i) - 2.0)) ;


    			t15 = lam(i) + psi(i) ;
    			t16 = psi(i)/t15 ;
    			t17 = pow(t16, psi(i)) ;
    			t18 = -psi(i) - 1.0 ;
    			t24 = lam(i)/t15 ;
    			t19 = -t24 + 1.0 ;
    			t20 = log(t19) ;
    			t21 = lam(i)*lam(i) ;
    			t22 = t15*t15 ;
    			t23 = -psi(i) - 2.0 ;
    			moments(i,3) = (t17*lam(i)*psi(i)*exp(t20*t18))/t15 - 3.0*t21*pow(t15, -2.0)*t17*t18 \
				* psi(i)*exp(t20*t23) + (t21*t23*t17*t18*lam(i)*psi(i)*exp(t20*(-psi(i) - 3.0)))/(t22*t15) ;

 	  			t16 = lam(i) + psi(i) ;
    			t17 = psi(i)/t16 ;
    			t18 = pow(t17, psi(i)) ;
    			t19 = -psi(i) - 1.0 ;
    			t25 = lam(i)/t16 ;
    			t20 = -t25 + 1.0 ;
    			t21 = log(t20) ;
    			t22 = lam(i)*lam(i) ;
    			t23 = t16*t16 ;
    			t24 = -psi(i) - 2.0 ;
    			t26 = -psi(i) - 3.0 ;
    			moments(i,4) = (t18*lam(i)*psi(i)*exp(t21*t19))/t16 - 7.0*t22*pow(t16, -2.0)*t18*t19 \
				*psi(i)*exp(t21*t24) + (6.0*t22*t24*t18*t19*lam(i)*psi(i)*exp(t21*t26))/(t23*t16) \
				-(t22*t22)*pow(t23, -2.0)*t24*t26*t18*t19*psi(i)*exp(t21*(-psi(i) - 4.0)) ;

    			t19 = lam(i) + psi(i) ;
    			t20 = psi(i)/t19 ;
    			t21 = pow(t20, psi(i)) ;
    			t22 = -psi(i) - 1.0 ;
    			t28 = lam(i)/t19 ;
    			t23 = -t28 + 1.0 ;
    			t24 = log(t23) ;
    			t25 = lam(i)*lam(i) ;
    			t26 = t19*t19 ;
    			t27 = -psi(i) - 2.0 ;
    			t29 = -psi(i) - 3.0 ;
    			t30 = t25*t25 ;
    			t31 = t26*t26 ;
    			t32 = -psi(i) - 4.0 ;
    			moments(i,5) = (t21*lam(i)*psi(i)*exp(t22*t24))/t19 - 15.0*t21*t22*t25*pow(t19, -2.0)\
					*psi(i)*exp(t24*t27) - 10.0*t21*t30*t22*pow(t26, -2.0)*t27*t29*psi(i) \
					*exp(t32*t24) + (25.0*t21*t22*t25*t27*lam(i)*psi(i)*exp(t24*t29))/(t26*t19) \
					+ (t21*t30*t22*t32*t27*t29*lam(i)*psi(i)*exp(t24*(-psi(i) - 5.0)))/(t31*t19) ;

    			t20 = lam(i) + psi(i) ;
    			t21 = psi(i)/t20 ;
    			t22 = pow(t21, psi(i)) ;
    			t23 = -psi(i) - 1.0 ;
    			t29 = lam(i)/t20 ;
    			t24 = -t29 + 1.0 ;
    			t25 = log(t24) ;
    			t26 = lam(i)*lam(i) ;
    			t27 = t20*t20 ;
    			t28 = -psi(i) - 2.0 ;
    			t30 = -psi(i) - 3.0 ;
    			t31 = t26*t26 ;
    			t32 = t27*t27 ;
    			t33 = -psi(i) - 4.0 ;
    			t34 = -psi(i) - 5.0 ;
    			moments(i,6) = (t22*lam(i)*psi(i)*exp(t23*t25))/t20 - 31.0*pow(t20, -2.0)*t22*t23*t26\
					*psi(i)*exp(t25*t28) - 65.0*t30*t22*t31*t23*pow(t27, -2.0)*t28*psi(i)*exp(t33 \
					*t25) + (90.0*t22*t23*t26*t28*lam(i)*psi(i)*exp(t30*t25))/(t20*t27) \
					+ (15.0*t30*t22*t31*t23*t33*t28*lam(i)*psi(i)*exp(t25*t34))/(t20*t32) \
					- (t30*t22*t31*t23*t33*t34*t26*t28*psi(i)*exp(t25*(-psi(i) - 6.0)))/(t32*t27) ;

    			t21 = lam(i) + psi(i) ;
    			t22 = psi(i)/t21 ;
    			t23 = pow(t22, psi(i)) ;
    			t24 = -psi(i) - 1.0 ;
    			t30 = lam(i)/t21 ;
    			t25 = -t30 + 1.0 ;
    			t26 = log(t25) ;
    			t27 = lam(i)*lam(i) ;
    			t28 = t21*t21 ;
    			t29 = -psi(i) - 2.0 ;
    			t31 = -psi(i) - 3.0 ;
    			t32 = t27*t27 ;
    			t33 = t28*t28 ;
    			t34 = -psi(i) - 4.0 ;
    			t35 = -psi(i) - 5.0 ;
    			t36 = -psi(i) - 6.0 ;
    			moments(i,7) = (t23*lam(i)*psi(i)*exp(t24*t26))/t21 - 63.0*pow(t21, -2.0)*t23*t24*t27 \
					*psi(i)*exp(t26*t29) - 350.0*t31*t23*t32*t24*pow(t28, -2.0)*t29*psi(i) \
					*exp(t34*t26) + (301.0*t23*t24*t27*t29*lam(i)*psi(i)*exp(t31*t26))/(t21*t28) \
					+ (140.0*t31*t23*t32*t24*t34*t29*lam(i)*psi(i)*exp(t26*t35))/(t21*t33) \
					- (21.0*t31*t23*t32*t24*t34*t35*t27*t29*psi(i)*exp(t26*t36))/(t33*t28) \
					+ (t31*t23*t32*t24*t34*t35*t27*t36*t29*lam(i)*psi(i)*exp(t26*(-psi(i) - 7.0)))/(t21*t33*t28) ;

				t22 = lam(i) + psi(i) ;
    			t23 = psi(i)/t22 ;
    			t24 = pow(t23, psi(i)) ;
    			t25 = -psi(i) - 1.0 ;
    			t31 = lam(i)/t22 ;
    			t26 = -t31 + 1.0 ;
    			t27 = log(t26) ;
    			t28 = lam(i)*lam(i) ;
    			t29 = t22*t22 ;
    			t30 = -psi(i) - 2.0 ;
    			t32 = -psi(i) - 3.0 ;
    			t33 = t28*t28 ;
    			t34 = t29*t29 ;
    			t35 = -psi(i) - 4.0 ;
    			t36 = -psi(i) - 5.0 ;
    			t37 = -psi(i) - 6.0 ;
    			t38 = -psi(i) - 7.0 ;
    			moments(i,8) = (t24*lam(i)*psi(i)*exp(t25*t27))/t22 - 127.0*pow(t22, -2.0)*t24*t25*t28 \
					*psi(i)*exp(t30*t27) - 1701.0*t30*t32*t24*t33*t25*pow(t29, -2.0)*psi(i)*exp(t35\
					*t27) + (966.0*t30*t24*t25*t28*lam(i)*psi(i)*exp(t32*t27))/(t22*t29) + (1050.0*t30 \
					*t32*t24*t33*t25*t35*lam(i)*psi(i)*exp(t27*t36))/(t22*t34) - (266.0*t30*t32*t24 \
					*t33*t25*t35*t36*t28*psi(i)*exp(t27*t37))/(t34*t29) + (28.0*t30*t32*t24*t33 \
					*t25*t35*t36*t28*t37*lam(i)*psi(i)*exp(t27*t38))/(t22*t34*t29) - t30*t32*t24*(t33\
					*t33)*t25*pow(t34, -2.0)*t35*t36*t37*t38*psi(i)*exp(t27*(-psi(i) - 8.0)) ;

    			t25 = lam(i) + psi(i) ;
    			t26 = psi(i)/t25 ;
    			t27 = pow(t26, psi(i)) ;
    			t28 = -psi(i) - 1.0 ;
    			t34 = lam(i)/t25 ;
    			t29 = -t34 + 1.0 ;
    			t30 = log(t29) ;
    			t31 = lam(i)*lam(i) ;
    			t32 = t25*t25 ;
    			t33 = -psi(i) - 2.0 ;
    			t35 = -psi(i) - 3.0 ;
    			t36 = t31*t31 ;
    			t37 = t32*t32 ;
    			t38 = -psi(i) - 4.0 ;
    			t39 = -psi(i) - 5.0 ;
    			t40 = -psi(i) - 6.0 ;
    			t41 = -psi(i) - 7.0 ;
    			t42 = t36*t36 ;
    			t43 = t37*t37 ;
    			t44 = -psi(i) - 8.0 ;
    			moments(i,9) = (t27*lam(i)*psi(i)*exp(t30*t28))/t25 - 255.0*t31*pow(t25, -2.0)*t27*t28 \
					*psi(i)*exp(t30*t33) - 7770.0*pow(t32, -2.0)*t33*t35*t27*t36*t28*psi(i)*exp(t30\
					*t38) - 36.0*t40*t41*t33*t42*t35*t27*t28*pow(t37, -2.0)*t38*t39*psi(i)*exp(t30 \
					*t44) + (3025.0*t31*t33*t27*t28*lam(i)*psi(i)*exp(t30*t35))/(t32*t25) + (6951.0\
					*t33*t35*t27*t36*t28*t38*lam(i)*psi(i)*exp(t30*t39))/(t25*t37) - (2646.0*t31*t33\
					*t35*t27*t36*t28*t38*t39*psi(i)*exp(t30*t40))/(t32*t37) + (462.0*t31*t40*t33 \
					*t35*t27*t36*t28*t38*t39*lam(i)*psi(i)*exp(t30*t41))/(t32*t25*t37) + (t40*t41*t33\
					*t42*t35*t44*t27*t28*t38*t39*lam(i)*psi(i)*exp(t30*(-psi(i) - 9.0)))/(t25*t43) ;
		}
	}



	// normalization factor and shaping polynomial
	normfactor = normfactor.fill(0.0);
	poly = poly.fill(0.0);
	prediction = prediction.fill(0.0);
	for(i = 0; i <= order; i++)
	{
		for(j = 0; j <= order; j++)
		{
			normfactor = normfactor + gam(i)*gam(j)*moments.column(i+j);
			prediction = prediction + gam(i)*gam(j)*moments.column(i+j+1);
		}

		for(j = 0; j < n; j++)
		{
			poly(j) = poly(j) + gam(i)*pow(y(j), (double) i);
		}
	}


	// log density and prediction
	for(i = 0; i < n; i++)
	{
		poly(i) = log(DBL_EPSILON + poly(i)*poly(i));
		t1 = DBL_EPSILON + psi(i)+lam(i);
		logdensity(i) = lgamma(y(i)+psi(i)) - lgamma(psi(i)) - lgamma(y(i)+1.0) \
			+ psi(i)*log(DBL_EPSILON + psi(i) / t1) + y(i)*log(DBL_EPSILON + lam(i) / t1);
		logdensity(i) = poly(i) + logdensity(i) - log(DBL_EPSILON+normfactor(i));
		prediction(i) = prediction(i) / normfactor(i);
	}

	f_return(0) = logdensity;
	f_return(1) = "na";
	f_return(2) = prediction;
	return octave_value_list(f_return);

}
