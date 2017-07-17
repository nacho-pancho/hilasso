#include <iostream>     // for debugging mostly
#include <cmath>       
#include <gsl/gsl_sf.h> // special functions
#include <gsl/gsl_cdf.h> // cumulative distribution functions
#include <gsl/gsl_integration.h> // cumulative distribution functions
#include "lookup.h"

#define MIN_P 1e-9
/*
 * CDF for the LG distribution given by the convolution of
 * a Gaussian r.v. with std. dev. sigma and a Laplacian with *inverse-scale* parameter
 * lambda, that is, L(x)=(\lambda/2) e^{-\lambda |x|}
 *
 * Taken from 
 *
 * G. Motta et al, "The iDUDE framework for 
 *
 */
double cdf_lg_P(double x, double lambda, double sigma) {
  double C1  = M_SQRT1_2/sigma; // 1/(sqrt(2)*sigma)
  double C2 = lambda*sigma*sigma;
  //#if DEBUG
  //std::cout << "C1=" << C1 << " C2="  << C2 << std::endl;
  //#endif
  return 0.5*( 1.0 + gsl_sf_erf(C1*x)
	       - 0.5*exp(lambda*(0.5*C2+x))*(gsl_sf_erf(C1*(x+C2))+1.0)
	       - 0.5*exp(lambda*(0.5*C2-x))*(gsl_sf_erf(C1*(x-C2))-1.0) );
}

/*
 * PDF for the LG distribution
 * Computed by me. The integral is easily derived (although quite lengthy)
 * by completing the squares.
 * Numerically unstable for lambda*sigma > 100
 */
double pdf_lg(double x, double lambda, double sigma) {
  double C1  = M_SQRT1_2/sigma; // 1/(sqrt(2)*sigma)
  double C2 = lambda*sigma*sigma;
  return (lambda/4.0)*
    (
     exp(  lambda*(x+0.5*C2) )*( 1 - gsl_sf_erf(C1*(x+C2)) ) +
     exp( -lambda*(x-0.5*C2) )*( 1 + gsl_sf_erf(C1*(x-C2)) ) 
     );
}

/*
 * PDF for the LG distribution/ GSL quadrature funtion interface.
 * Parameters are precomputed constants
 *
 * C1  = M_SQRT1_2/sigma
 * C2 = lambda*sigma*sigma
 * C3 = 0.5*C2
 * C4 = lambda/4
 */
double F_pdf_lg(double x, void* params) {
  double sigma  = ( (double*)params )[0];
  double C2  = ( (double*)params )[1];
  double C3  = ( (double*)params )[2];
  double lambda  = ( (double*)params )[3];
  if (fabs(x) >= 10.0*(0.5/lambda/lambda + sigma*sigma)) {
    return 0.0;
  } else {
    double C1 = M_SQRT1_2/sigma;
    return 0.25*lambda*
      (
       exp(  lambda*(x+C3) )*( 1 - gsl_sf_erf(C1*(x+C2)) ) +
       exp( -lambda*(x-C3) )*( 1 + gsl_sf_erf(C1*(x-C2)) ) 
       );
  }
}

parametric_lut* create_laplacian_parametric_lut(int maxx, double qx, int M) {
  double maxt = 2.0*sqrt(double(maxx)); // empirically, maximum value of the LG parameter I found was sqrt(x), this is a little more because it's cheap to compute, small storage, and it's safer
  double qt = 1.0/sqrt(double(M)); // universal two-step quantization step
  parametric_lut* lut = new parametric_lut(double(maxx),qx, maxt, qt);
  //int minx   = -maxx;
  int offset =  lut->get_offset();
  int tsize  = lut->get_tsize();
  //double theta = 0.0;
  double ** lut_data = lut->get_data();
  //
  // for theta=0, this is the delta at 0
  //
  sc_vector<double> v(lut_data[0],tsize);
  sc_fill(-log2(MIN_P),v);
  v[offset] = 0.0;
  //
  // for theta>-
  //
  double theta = qt;
  for (int i= 1 ; i < lut->get_ntables(); i++, theta += qt) {
     sc_vector<double> curr = sc_vector<double>(lut_data[i],tsize);
     int j=offset;
     double C; // normalization constant
     double xl;
     double xr=0.5*qx;
     C = curr[j] = (1.0-exp(-fabs(xr)/theta)); // x = 0 
     for (j=offset+1 ;j < lut->get_tsize(); j++) { // offset    
       xl =  xr;
       xr += qx;
       double P = (j < (lut->get_tsize()-1))? 
	 0.5*(exp(-fabs(xl)/theta) - exp(-fabs(xr)/theta)) :
	 0.5* exp(-fabs(xl)/theta); // xr = inf 
       curr[j] = P;
       C += 2.0*P; // 0 is counted once, the rest twice
       if (P < MIN_P)
 	break;
     }
     // tail is zero
     for (; j < lut->get_tsize(); j++)  {
       // numerically, I dont want Inf
       curr[j] = MIN_P; 
       C += 2.0*MIN_P;
     }
     // 
     // normalization and log
     //
     for (j=offset; j < tsize; j++) {
       curr[2*offset-j] = curr[j] = log2(C)-log2(curr[j]);
     }
   }
   return lut;
}

parametric_lut* create_moe_parametric_lut(int maxx, double qx, int M, double kappa) {
  double maxt = 2.0*sqrt(double(maxx)); // empirically, maximum value of the LG parameter I found was sqrt(x), this is a little more because it's cheap to compute, small storage, and it's safer
  double qt = 1.0/sqrt(double(M)); // universal two-step quantization step
  parametric_lut* lut = new parametric_lut(double(maxx),qx, maxt, qt);
  //int minx   = -maxx;
  int offset =  lut->get_offset();
  int tsize  = lut->get_tsize();
  //double theta = 0.0;
  double ** lut_data = lut->get_data();
  //
  // for theta=0, this is the delta at 0
  //
  sc_vector<double> v(lut_data[0],tsize);
  sc_fill(-log2(MIN_P),v);
  v[offset] = 0.0;

  //
  // for theta>-
  //
  double theta = qt;
  for (int i= 1 ; i < lut->get_ntables(); i++, theta += qt) {
     double beta = (kappa-1.0)*theta;
     sc_vector<double> curr = sc_vector<double>(lut_data[i],tsize);
     int j=offset;
     double C; // normalization constant
     double xl;
     double xr=0.5*qx;
     C = curr[j] = (1.0 - pow(fabs(xr)/beta+1.0 , -kappa)); // x = 0 
     for (j=offset+1 ;j < lut->get_tsize(); j++) { // offset    
       xl =  xr;
       xr += qx;
       double P = (j < (lut->get_tsize()-1))? 
	 pow(fabs(xl)/beta+1.0,-kappa) - pow(fabs(xr)/beta+1.0,-kappa) :
	 pow(fabs(xl)/beta+1.0,-kappa); // xr = inf
       curr[j] = P;
       C += 2.0*P; // 0 is counted once, the rest twice
       if (P < MIN_P)
 	break;
     }
     // tail is zero
     for (; j < lut->get_tsize(); j++)  {
       // numerically, I dont want Inf
       curr[j] = MIN_P; 
       C += 2.0*MIN_P;
     }
     // 
     // normalization and log
     //
     for (j=offset; j < tsize; j++) {
       curr[2*offset-j] = curr[j] = log2(C)-log2(curr[j]);
     }
   }
   return lut;
}

parametric_lut* create_lg_parametric_lut(int maxx, double qx, int M, double sigma) {
  double maxt = 2.0*sqrt(double(maxx)); // empirically, maximum value of the LG parameter I found was sqrt(x), this is a little more because it's cheap to compute, small storage, and it's safer
  double qt = 1.0/sqrt(double(M)); // universal two-step quantization step
  parametric_lut* lut = new parametric_lut(double(maxx),qx, maxt, qt);
  int minx   = -maxx;
  int offset =  maxx;
  int tsize  = lut->get_tsize();
  // we use quadrature functions from GSL
  //
  //
  gsl_function F;
  double params[4];
  F.function = &F_pdf_lg;
  F.params   = &params[0];
  // 0: C1(sigma)  = M_SQRT1_2/sigma
  // 1: C2(sigma,lambda) = lambda*sigma*sigma = sigma*sigma/theta
  // 2: C3(sigma,lambda) = 0.5*C2
  // 3: C4(lambda) = lambda

  //
  // compute lookup function for each theta. For theta=0 this is only
  // a discretized Gaussian.
  // for theta < certain value, we also use a Gaussian because the
  // LG functions are at their current state very unstable for very small theta
  // 0.25*20 <= 5 
  // sigma/theta <= 5 -> theta >= sigma/5
  //
  params[0] = sigma;
  double theta = 0.0;
  double min_theta = sigma/5.0;
  int i;
  double** lut_data = lut->get_data();
  for (i = 0; theta < min_theta; i++, theta += qt) {
    int x;
    double C = 0.0;
    sc_vector<double> curr = sc_vector<double>(lut_data[i],tsize);
    for (x=0; x <= maxx; x++) { // offset    
      double Fl = (x>minx) ? gsl_cdf_gaussian_P(double(x)-0.5,sigma): 0.0;
      double Fr = (x<maxx) ? gsl_cdf_gaussian_P(double(x)+0.5,sigma): 1.0;
      double P = Fr-Fl;
      curr[offset+x] = P;
      C += (x>0)? 2.0*P: P; // 0 is counted once, the rest twice
      if (P < MIN_P)
	break;
    }
    //
    // fill up tails with negligible nonzero probability
    //
    for (; x<= maxx; x++)  {
      // numerically, I dont want Inf
      curr[offset+x] = MIN_P; 
      C += 2.0*MIN_P;
    }
    //
    // normalize
    //
    for (x=0; x<=maxx; x++) {
      curr[offset-x] = curr[offset+x] = log2(C)-log2(curr[offset+x]);
    }
  }
  // Note that here theta is the scale parameter of the Laplacian,
  // that is, L(x;theta) = \frac{1}{2\theta} e^{-\frac{|x|}{\theta}}
  // but the LG parameter that we pass is \lambda=\frac{1}{\theta}
  // 
  gsl_integration_workspace* wks = gsl_integration_workspace_alloc(1000);
  for ( ; i < lut->get_ntables(); i++, theta += qt) {
    
    sc_vector<double> curr = sc_vector<double>(lut_data[i],tsize);
    
    params[1] = sigma*sigma/theta;
    params[2] = 0.5*params[1];
    params[3] = 1.0/theta;
    int x;
    double C=0.0;
    for (x=0; x <= maxx; x++) { // offset    
      double xl = (x>minx) ? double(x)-0.5 : -1e4;
      double xr = (x<maxx) ? double(x)+0.5 :  1e4;
      double abserr;
      double P;
      gsl_integration_qag(&F,        // pointer to gsl_function
			  xl, // integration interval [a,...
			  xr, //  b]
			  1e6, // abs error
			  1e-3, // rel error
			  1000, // maximum number of intervals for computation
			  GSL_INTEG_GAUSS15, // number of samples in Gauss-Kronrod technique: 15
			  wks,
			  &P, &abserr);
      curr[offset+x] = P;
      C += (x>0)? 2.0*P: P;
      if (P < MIN_P)
	break;
    }
    //
    // fill up tails with negligible nonzero probability
    //
    for (; x<= maxx; x++)  {
      // numerically, I dont want Inf
      curr[offset+x] = MIN_P;
      C += 2.0*MIN_P; 
    }
    //
    // normalize and take log
    //
    for (x=0; x<=maxx; x++) {
      curr[offset-x] = curr[offset+x] = log2(C)-log2(curr[offset+x]);
    }
  } // end for each theta
  gsl_integration_workspace_free(wks);  
  
  return lut; 
}

/**
 * by exchanging integration order, each table here is simply a mixture
 * of the above, weighted by the mixing prior of MOE w(theta;beta,kappa)
 * This we do by creating an auxiliary LG LUT and then a Gamma weights vector for each
 * parameter, and then just do matrix/vector multiplication!
 */
parametric_lut* create_moeg_parametric_lut(int maxx, double qx, int M, double kappa, double sigma) {
  parametric_lut *lg_lut, *moeg_lut;
  double *w_data;
  //
  // auxiliary LG lut: we will average the LG numerically using this
  //
  double qt = 1.0/sqrt(double(M)); // universal two-step quantization step
  double maxt = 2.0*sqrt(double(maxx)); // empirically, maximum value of the LG parameter I found was s

  lg_lut = create_lg_parametric_lut(maxx,qx,M,sigma);
  moeg_lut = new parametric_lut(maxx,qx,maxt,qt);

  w_data = new double[lg_lut->get_ntables()];
  sc_vector<double> w(w_data,lg_lut->get_ntables());
  //
  // beta = 0, the weighting function is a delta at 0
  // this means that the first 'mixture' is a Gaussian, 
  // same as the first one in LG
  //
  double** lg_lut_data = lg_lut->get_data();
  double** moeg_lut_data = moeg_lut->get_data();
  sc_vector<double> u(lg_lut_data[0],lg_lut->get_tsize());
  sc_vector<double> v(moeg_lut_data[0],lg_lut->get_tsize());
  sc_copy(u,v);
  //
  // theta > 0
  //   
  //
  // undo the -log2 so that we get back probability tables
  //
  for (int i = 0; i < lg_lut->get_ntables(); i++) {
    for (int j = 0; j < lg_lut->get_tsize(); j++) {
      lg_lut_data[i][j] = pow(2.0,-lg_lut_data[i][j]);
    }
  }

  for (int i = 1; i < moeg_lut->get_ntables(); i++) {    
    double beta = (kappa-1.0)*qt*i;
    //
    // instantiate gamma 
    //
    double Fl(0.0),Fr(0.0);
    double tl = 0.0;
    double tr = lg_lut->get_qt()*0.5; 
    int j;
    for (j=0; j< (lg_lut->get_ntables()-1) ; j++) {    
      Fl = gsl_cdf_gamma_P(tl,kappa,beta);
      Fr = gsl_cdf_gamma_P(tr,kappa,beta);
      w[j] = Fr-Fl;
      tl = tr;
      tr += lg_lut->get_qt();
    }
    // last one
    w[j] = 1.0 - Fr;
    //std::cout << "sum(w)=" << sc_sum(w) << std::endl;
    //
    // w(.;beta) * lg_lug = moeg_lut(.;beta)
    //
    sc_vector<double> Q(moeg_lut_data[i],moeg_lut->get_tsize());
    sc_fill(0.0,Q);
    for (j=0; j < lg_lut->get_ntables(); j++) {
      sc_vector<double> P(lg_lut_data[j],lg_lut->get_tsize());
      sc_scaled_add(w[j],P,Q);
    }
    //
    // take -log2 again
    //
    for(j=0; j < moeg_lut->get_tsize(); j++)
      moeg_lut_data[i][j] = -log2(moeg_lut_data[i][j]);
  }
  

  delete[] w_data;
  delete lg_lut;
  return moeg_lut;
}

