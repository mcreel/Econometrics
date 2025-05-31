
#include <oxstd.h>
#import <maximize>

static decl y;   		//	count dep. vbl.
static decl x;			//  used to get starting values in negbin and logit
static decl add_up;		//  sum obj. or not
static decl negbin_type;


#include "NegBinSNP.ox"
#include "MLE.ox"


main()
{
	decl data, n, model, theta, names, title, obj_value, p;
	
	p = 2; // degree of polynomial
// Read in the data
	data = loadmat("health.mat");

// Define dep and expl vbls
	/* 
	The dep. vbls, in corresponding 	column
	OBDV: Office based doctor visits	    0
	OPV: Outpatient doctor visits			1
	ERV: Emergency room visits				2
	IPV: Inpatient visits					3
	DV: Dental visits						4
	PRESCR: Prescriptions					5
	*/

// all 6 dependent variables	
	title = "MEPS data, OBDV";
	y = data[][0];

	/*
	The available expl. variables and their columns
    PUBLIC_INS  	6
	PRIVATE_INS		7
	SEX				8
	AGE				9
	EDUC			10
	INCOME			11
	*/

	x = data[][6:11];
	x = standardize(x);
	n = rows(x);
	x = ones(n,1) ~ x;

	names = {"constant","pub_ins","priv_ins","sex","age","educ","inc","ln_alpha"};
	model = "negbin_snp_obj";
	negbin_type = 1;
	theta = zeros(sizec(x)+1,1);
	theta = theta | zeros(p,1);
	theta = mle(model, theta, names, title, 1);
	
}
