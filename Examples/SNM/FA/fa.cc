// factor ARCH: Billio and Monfort's model

#include <oct.h>
#include <octave/Cell.h>

DEFUN_DLD(fa, args, ,"fa")
{
	// parameter of model
	ColumnVector model_params (args(0).column_vector_value());
	double a1 = model_params(0);
	double a2 = model_params(1);
	double sig = model_params(2);
	double b2 = model_params(3);


	// other args
	Cell modelargs (args(1).cell_value());
	Matrix randdraws (modelargs(0).matrix_value()); // random numbers, passed this way to keep fixed
	int n (modelargs(1).int_value());
	int burnin (modelargs(2).int_value()); // number of burnin periods to eliminate startup effects

	int t;
	double y1, y2, ystar, ystarlag, sig_ystar;

	Matrix ys(n,2);
	ys.fill(0.0);
	octave_value_list f_return;
	ystar = 0.0;
	ystarlag = 0.0;
	// main loop, first burnin, then keepers
	for (t=0; t < n + burnin; t++) { // loop over time
		sig_ystar = sqrt(a1 + a2*ystarlag*ystarlag);
		ystar = sig_ystar * randdraws(t,0);
		y1 = ystar + sig * randdraws(t,1);
		y2 = b2*ystar + sig * randdraws(t,2);
		if (t >= burnin) {
			ys(t-burnin,0) =  y1;
			ys(t-burnin,1) =  y2;
		}
		ystarlag = ystar;
	}

	f_return(0) = ys;
	return f_return;
}
