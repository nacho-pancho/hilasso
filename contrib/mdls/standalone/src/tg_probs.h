#ifndef _TG_PROBS_
#define _TG_PROBS_
#include <gsl/gsl_cdf.h>
#include <math.h>
//#include "Counts.h"


double discrete_gaussian_P(int x, double sigma);
double TSGD_P(int x, double theta, double d);
#ifdef SLOW
double TSGD_Gauss_P(int z, double theta, double d, double sigma);
#else
#define TSGD_Gauss_P tsgd_gauss_discrete
#endif
/**
 * evaluate cumulative TSGD+GAUSS distribution
 */
double tsgd_gauss_cumul(double z, double theta, double d, double sigma);
/**
 * evaluate a single value from the discrete T+G
 */
double tsgd_gauss_discrete(int z, double theta, double d, double sigma);
/**
 * compute the full T+G for a given support [min,max] into a vector pmf
 */
//void tsgd_gauss_pmf(double theta, double d, double sigma,info::Counts& pmf, double min_prob, bool accu_tails);

double tsgd_gauss_density(double z, double theta, double d, double sigma);
#endif

