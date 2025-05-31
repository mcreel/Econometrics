// logarithmic stochastic volatility model
// used by Andersen et al., Chumacero, Takeda, etc.
// y(t) = exp(ystar(t)/2) * e1
// ystar(t) = a + b  ystar(t-1) + sig*e2
//
// a different parameterization (e.g., Fermanian and Salanie, Altissimo and Mele) is
// y(t) = exp(a/2)*exp(ystar(t)/2) * e1
// ystar(t) = b*ystar(t-1) + sig*e2
//
// the difference is whether a or exp(a) is estimated. Other than that, they are equivalent

#include <oct.h>
#include <octave/Cell.h>

DEFUN_DLD(sv, args, ,"sv")
{
	// parameter of model
	ColumnVector model_params (args(0).column_vector_value());
	double phi = model_params(0);
	double sigb = model_params(1);
	double sige = model_params(2);

	// other args
	Cell modelargs (args(1).cell_value());
	Matrix randdraws (modelargs(0).matrix_value()); // random numbers, passed this way to keep fixed
	int n (modelargs(1).int_value());
	int burnin (modelargs(2).int_value()); // number of burnin periods to eliminate startup effects

	int t;
	double y, ystar;
	double ystarlag = 0.0;

	ColumnVector ys(n);
	ys.fill(0.0);
	octave_value_list f_return;

	// main loop, first burnin, then keepers
	for (t=0; t < (n + burnin); t++) { // loop over time
		ystar = phi*ystarlag + sige*randdraws(t,1);
		y = sigb*exp(ystar / 2.0) * randdraws(t, 0);
		if (t >= burnin) {
			ys(t-burnin) +=  y;
		}
		ystarlag = ystar;
	}
	f_return(0) = ys;
	return f_return;
}
