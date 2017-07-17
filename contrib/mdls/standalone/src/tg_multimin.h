#ifndef _MULTIMIN_
#define _MULTIMIN_

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_vector.h>
#include "tsgd+gauss.h"

#define MULTIMIN_N  3 // number of variables in the minimization
#define MULTIMIN_SIZE_TARGET 5e-4  // threshold for stopping minimization
#define MULTIMIN_MAX_ITERS   150   // MAX No. of iterations allowed

class multimin_params {
public:
	iterbounds t;  // bounds for theta
	iterbounds dd; // bounds for d
	iterbounds s;  // bounds for sigma
	cost_type cost;
	histogram *hp;

	multimin_params() {}
	multimin_params(cost_type c, iterbounds tt, iterbounds ddd, iterbounds ss, histogram& h)
		:t(tt),dd(ddd),s(ss),cost(c),hp(&h) {}


};

double tg_minimize_cost_solver(cost_type cost, histogram& h, 
				  const tdsvec& init_point, 
				  tdsvec& result_point, 
				  iterbounds t, iterbounds dd, iterbounds s);

double vartsgd(double theta, double d);
double thetafromvar(double var, double d);

#endif
