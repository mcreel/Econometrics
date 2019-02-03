 
#include <oct.h>
#include <octave/parse.h>
#include <octave/Cell.h>
#include <octave/lo-mappers.h>
#include <float.h>
#include "error.h"

DEFUN_DLD(emm_moments_C, args, , "emm_moments, C++ version")
{
  int nargin = args.length();
  if (nargin != 3)
  {
    error("emm_moments: you must provide 3 arguments");
    return octave_value_list();
  }

//   // check the arguments
//   if (any_bad_argument(args)) return octave_value_list();

	Matrix theta (args(0).matrix_value());
	Matrix data (args(1).matrix_value());
  Cell momentargs (args(2).cell_value());

	// the things in momentargs
	int k (momentargs(0).int_value());
  std::string dgp (momentargs(1).string_value());
 	Cell dgpargs (momentargs(2).cell_value());
  std::string sg (momentargs(3).string_value());
  Cell sgargs (momentargs(4).cell_value());
	Matrix phi (momentargs(5).matrix_value());

	int g = phi.rows(); // number of moments
	int n = data.rows();  // number of observations
	int c = data.columns(); // columns of data
	Matrix sgdata = data.extract_n(0, 0, n-1, k); // y and x together (we'll replace y with simulated)
	Matrix x = data.extract_n(0, 1, n-1, k); // next cols are exogs

	// random draws passed in data to ensure fixed over estimation
	Matrix rand_draws = data.extract_n(0, k+1, n-1, c-1); // remaining cols are random draws

	Matrix scores(n, g, 0.0); // container for moment contribs
// 	int reps = rand_draws.columns(); // number of simulations
// 	Matrix y, e; 
// 	int i;
  	octave_value_list f_return; // holder for feval returns
// 
//   octave_value_list g_args(3,1); // to evaluate numgradient using celleval
//   g_args(0) = sg;
//   g_args(2) = sgargs;
// 	
//   octave_value_list d_args(5,1); // to evaluate dgp using celleval
//   d_args(0) = dgp;
// 	d_args(1) = theta;
// 	d_args(2) = x;
// 	d_args(4) = dgpargs;
// 	
// 	
// 	for(i = 0; i < reps; i++)
// 	{
// 		e = rand_draws.extract_n(0, i, n-1, 1);
// 		d_args(3) = e;
//     f_return = feval("celleval", d_args);
// 		y = f_return(0).matrix_value();
// 		sgdata.insert(y,0,0);
// 		g_args(1) = sgdata;
// 		f_return = feval("numgradient", g_args);
//  		scores = scores + f_return(0).matrix_value();
//  	}
//  		//	scores = scores /// reps; // average over number of simulations

	f_return(0) = scores;
	return f_return;
}
