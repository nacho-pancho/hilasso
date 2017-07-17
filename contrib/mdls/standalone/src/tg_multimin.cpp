#include <stdio.h>

#include <iostream>
#include <iomanip>
#include <time.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include "tsgd+gauss.h"
#include "tg_multimin.h"

using namespace std;
gsl_multimin_function my_func;

double f_cost(const gsl_vector *tds, void *params)
{
	double t = gsl_vector_get(tds,0), 
		   d = gsl_vector_get(tds,1), 
		   s = gsl_vector_get(tds,2);
	histogram *hp = ((multimin_params *)params)->hp;

	cost_type cost = ((multimin_params *)params)->cost;
	if ( t < MIN_THETA ) t = MIN_THETA; 
	else if ( t > MAX_THETA ) t = MAX_THETA;  // don't take invalid values of theta
	double c = hp->cost_fxn(cost, t, d,  s);
	return c;
}

void df_cost(const gsl_vector *tds, void *params, gsl_vector *g) {
   double t = gsl_vector_get(tds,0), 
		  d = gsl_vector_get(tds,1), 
		  s = gsl_vector_get(tds,2);
	histogram *hp = ((multimin_params *)params)->hp;

	cost_type cost = ((multimin_params *)params)->cost;

	tdsvec gr;
	hp->tg_gradient(cost, t, d, s, gr);
	gsl_vector_set(g, 0, gr[0]);
	gsl_vector_set(g, 1, gr[1]);
	gsl_vector_set(g, 2, gr[2]);
}


void fdf_cost (const gsl_vector *x, void *params,
             double *f, gsl_vector *df)
{
       *f = f_cost(x, params);
       df_cost(x, params, df);
}

void normalize_gsl_tds(gsl_vector *x, multimin_params params)
{
	double xt, xd, xs;
	xt = gsl_vector_get(x, 0);
	xd = gsl_vector_get(x, 1);
	xs = gsl_vector_get(x, 2);

	params.t.enforce(xt);
	params.dd.enforce(xd);
	params.s.enforce(xs);

	gsl_vector_set(x, 0, xt);
	gsl_vector_set(x, 1, xd);
	gsl_vector_set(x, 2, xs);
}



double tg_minimize_cost_solver(cost_type cost, histogram &h, 
				  const tdsvec& init_point, tdsvec& result_point, 
				  iterbounds t, iterbounds dd, iterbounds s)  // bounds on variables
{
	
	multimin_params params(cost, t, dd, s, h);  // build parameter structure

	const gsl_multimin_fminimizer_type *T;
	gsl_multimin_fminimizer *fm;
    
	// Starting point
	gsl_vector *x = gsl_vector_alloc(MULTIMIN_N);

	init_point.to_gsl_vec(x);
	     
    my_func.f = &f_cost;
    my_func.n = MULTIMIN_N;
    my_func.params = (void *)(&params);	

	T = gsl_multimin_fminimizer_nmsimplex;
    fm = gsl_multimin_fminimizer_alloc (T, MULTIMIN_N);

    
	//initialize minimizer
    //	const double tolerance = 0.01;
    //	const double step_size = 1e-4;
	gsl_vector *v_step_size = gsl_vector_alloc(MULTIMIN_N);
	gsl_vector_set(v_step_size, 0, t.i);
	gsl_vector_set(v_step_size, 1, dd.i);
	gsl_vector_set(v_step_size, 2, s.i);
	gsl_multimin_fminimizer_set (fm, &my_func, x,  v_step_size);

	int iter = 0;
	int status;
	double fval;
	do {
          iter++;

		  // normalize input		  
          status = gsl_multimin_fminimizer_iterate (fm);
		  normalize_gsl_tds(fm->x, params);  // normalize x within boundaries
     
          if (status)
             break;
		   status = gsl_multimin_test_size (fm->size, MULTIMIN_SIZE_TARGET);
#if 0    
          if (status == GSL_SUCCESS)  printf ("Minimum found at:\n");
		  fval = gsl_multimin_fminimizer_minimum(fm);
     
          printf ("%5d %.5f %.5f %.5f %10.5f\n", iter,
                   gsl_vector_get (fm->x, 0),
                   gsl_vector_get (fm->x, 1),
		  	     gsl_vector_get (fm->x, 2),
                   fval);
#endif     
      } while (status == GSL_CONTINUE && iter < MULTIMIN_MAX_ITERS);
	  cout << "Multimin iterations: " << iter << endl;
	  result_point.from_gsl_vec(fm->x);
      fval = gsl_multimin_fminimizer_minimum(fm);
      gsl_multimin_fminimizer_free (fm);
      gsl_vector_free (x);
     
      return fval;
}

double vartsgd(double theta, double d) {
	double cg = theta * (4 * pow(theta, (double) (2 * d + 1)) 
		+ theta * theta + pow(theta, (double) (2 + 2 * d)) + 
		pow(theta, (double) (2 * d)) + pow(theta, (double) (4 * d))) 
		* pow(pow(theta, (double) (2 * d)) + theta, -2) * pow(theta - 1, -2);
	return cg;
}

struct tsgd_param {
	double var;
	double d;
};

double my_vartsgd(double theta, void *p) {
	struct tsgd_param *param = (struct tsgd_param *)p;
	double d = param->d;
	double var = param->var;
	return (vartsgd(theta, d)-var);
}

#define _TG_SOLVER_TOLERANCE 1e-6
double thetafromvar(double var, double d)
{
	try {
		const gsl_root_fsolver_type * T  = gsl_root_fsolver_bisection;
		gsl_root_fsolver * s   = gsl_root_fsolver_alloc (T);
		gsl_function F;
		F.function = &my_vartsgd;
		int iter = 0, status, end_cond = 0;

		struct tsgd_param param;
		param.d = d;
		param.var = var;

		F.params = &param;
		double x_lo = 0.01, x_hi = 0.99;

		double y_lo = my_vartsgd(x_lo, F.params);
		double y_hi = my_vartsgd(x_hi, F.params);

		if ( !(y_lo*y_hi < 0) ) {
			cerr << "*** Bad root solver initialization" << endl;
			throw 10;
		}

		//int rc = gsl_root_fsolver_set(s, &F, x_lo, x_hi);
		gsl_root_fsolver_set(s, &F, x_lo, x_hi);

		do {
			iter++;
			status = gsl_root_fsolver_iterate(s);
			//double r = gsl_root_fsolver_root(s);
			gsl_root_fsolver_root(s);
			x_lo = gsl_root_fsolver_x_lower (s);
			x_hi = gsl_root_fsolver_x_upper (s);
			// cout << "iter " << iter << " lo: " << x_lo << " hi: " << x_hi << endl;
			end_cond = ( fabs(x_hi-x_lo) < _TG_SOLVER_TOLERANCE );
		} while ( !end_cond ) ;
		return (x_lo+x_hi)/2;
	}
	catch (int e) {
		cerr << "Root solver failed to find theta from variance. Error code: " << e << endl;
		return -1;
	}
}
