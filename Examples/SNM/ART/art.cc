// autoregressive Tobit - Fermanian and Salanié's model


#include <oct.h>
#include <octave/Cell.h>

DEFUN_DLD(art, args, ,"art")
{

	// parameter of model
	ColumnVector model_params (args(0).column_vector_value());
	double a = model_params(0);
	double b = model_params(1);
	double sig = model_params(2);

	// other args
	Cell modelargs (args(1).cell_value());
	Matrix randdraws (modelargs(0).matrix_value()); // random numbers, passed this way to keep fixed
	int n (modelargs(1).int_value());
	int burnin (modelargs(2).int_value()); // number of burnin periods to eliminate startup effects

	int t;
	double e, z, y;
	double zlag = 0.0;

	ColumnVector ys(n);
	ys.fill(0.0);
	octave_value_list f_return;

	// main loop, first burnin, then keepers
	for (t=0; t < (n + burnin); t++) { // loop over time
		e = sig * randdraws(t,0);
		// latent variable
		z = a + b*zlag + e;
		// observed variable
		y = 0.0;
		if (z > 0.0) {
			y = z;
		}
		if (t >= burnin) {
			ys(t-burnin) +=  y;
		}
		zlag = z;
	}

	f_return(0) = ys;
	return f_return;
}
