// Copyright (C) 2005,2007  Michael Creel <michael.creel@uab.es>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// NOTE: in this program, k is lagged k. That is k(t) is the value of capital in period t-1
// this is just to simplify the program
#include <oct.h>
#include <octave/Cell.h>
#include <octave/oct-rand.h>

DEFUN_DLD(simple_growth, args, ,"simple_growth.cc: example model for solution by pea.m\n\
This is illustrative of the i/o format the model needs to have.\n\
This simulates the simple stochastic growth model.")
{
	// extract the arguments
	ColumnVector exp_params (args(0).column_vector_value());
	ColumnVector model_params (args(1).column_vector_value());
	double randseed (args(2).double_value());
	Cell pea_args (args(3).cell_value());
	double up_bound (args(4).double_value());
	double low_bound (args(5).double_value());

	// parameters of model
	double alpha = model_params(0);
	double delta = model_params(1);
	double gam = model_params(2);
	double d = model_params(3);
	double sigma = model_params(4);
	double rho = model_params(5);

	// PEA controls
	int T = (pea_args(0).int_value());
	int burnin = (pea_args(1).int_value());
	int nodes = (pea_args(5).int_value());
	nodes = nodes + 1;

	// set up RNG and get draws
	octave_rand::seed(randseed); // set the seed to the fixed value to eliminate chatter
	octave_rand::distribution("normal");  // we'll be using draws from U(0,1)

	// other declares and initializations
	int t, tt, hit;
	ColumnVector k(burnin + T + 1);
	ColumnVector c(burnin + T);
	ColumnVector shocks(burnin + T);
	Matrix data(T-1,4);;
	double uprime, epsi;
	octave_value_list f_return; // holder for returns (k, e, hit)
	data.fill(1.0); // second column is constant

	// use steady state to start up
	k(0) = pow(((1-delta+delta*d) / (alpha*delta)) , (1/(alpha-1)));
	hit = 0;

	// simulate, with burnin
	for (t=0; t < (burnin+T); t++)
	{
		if (t == 0) {
			shocks(t) = 1.0;
		}
		else {
			epsi = octave_rand::scalar();
			epsi = sigma*epsi;
			shocks(t) = pow(shocks(t-1), rho) * exp(epsi);
		}
		uprime = exp(exp_params(0) + exp_params(1)*log(k(t)) + exp_params(2)*log(shocks(t)));
		c(t) = pow(delta*uprime, (-1/gam));
		k(t+1) = pow(k(t), alpha) * shocks(t) - c(t) + (1-d)*k(t);
		if (k(t+1) > up_bound) {
			k(t+1) = k(0); hit = 1;
		}
		if (k(t+1) < low_bound) {
    			k(t+1) = k(0); hit = 1;
		}
		c(t) = pow(k(t), alpha) * shocks(t)  + (1-d)*k(t)-k(t+1);

	}

	// chop off the burnin period
	k = k.extract(burnin,k.rows());
 	c = c.extract(burnin,c.rows());
 	shocks = shocks.extract(burnin, shocks.rows());

	// data for model of expectations
	for (t = 0 ; t < data.rows(); t++) {
		data(t,0) = pow(c(t+1),(-gam)) * (1 - d + pow(k(t+1),(alpha-1)) * alpha * shocks(t+1));
		data(t,2) = log(k(t));
		data(t,3) = log(shocks(t));
	}
	f_return(1) = hit;
	f_return(0) = data;
  	return f_return;
}
