// dgp for AR1 model

#include <oct.h>
#include <octave/Cell.h>

DEFUN_DLD(ar1, args, ,"ar1")
{
	// arguments are rho, N, T

	ColumnVector model_params (args(0).column_vector_value());
	double rho0 = model_params(0);
	double rho1 = model_params(1);

	// other args
	Cell modelargs (args(1).cell_value());
	Matrix randdraws (modelargs(0).matrix_value()); // random numbers, passed this way to keep fixed
	int n (modelargs(1).int_value());
	int burnin (modelargs(2).int_value()); // number of burnin periods to eliminate startup effects


	// declares
	int t;
	ColumnVector ys(n);
	ys.fill(0.0);
	double y, ylag;
	ylag = 0.0;
	octave_value_list f_return;

	// main loop, first burnin, then keepers
	for (t=0; t < (n + burnin); t++) { // loop over time
		y = rho0 + rho1*ylag + randdraws(t, 0);
		if (t >= burnin) {
			ys(t-burnin) +=  y;
		}
		ylag = y;
	}
	f_return(0) = ys;
	return f_return;
}
