#include "thresholding.h"

//
// returns true if converged, false if max_iter reached
//
template<class T>
bool hilasso_admom(const T* x, 
		   int n, 
		   const T& lambda1,
		   const T& lambda2,
		   int max_iter,
		   const T& tol,
		   const T& c,
		   T* b)
{
  T p[n], beta[n], xhat[n]; // lagrange multiplier, aux variable
  T c1,c2;
  c1 = T(1)/(c+T(1));
  c2 = T(1)/(c);
  //  
  // initialization
  //
  //printf("n=%d\tlambda1=%g\tlambda2=%g\tmax_iter=%d\ttol=%g\tc=%g\n",
  //	  n,lambda1,lambda2,max_iter,tol,c);
  for (int i=0; i < n; i++) {
    beta[i] = b[i] = T(0);
    p[i] = T(1);
  }
  //
  // main loop
  //
  T dif;
  int j;
  for (j=0; j < max_iter; j++) {
    //
    // update b
    //
    // S(x+c*beta-p,lambda1)
    for (int i =0; i < n; i++) {
      xhat[i] = x[i]+c*beta[i]-p[i];
    }
    scalar_soft_thres(xhat,n,lambda1,b);
    //
    // b = b*(1/(c+1)):
    //
    for (int i =0; i < n; i++) {
      b[i] = c1*b[i];
    }    
    //
    // update beta
    //
    // S(p+c*b,lambda2)
    //
    for (int i =0; i < n; i++) {
      xhat[i] = p[i]+c*b[i];
    }
    vector_soft_thres(xhat,n,lambda2,beta);
    //
    // beta = beta*(1/c):
    //
    for (int i =0; i < n; i++) {
      beta[i] = c2*beta[i];
    }    
    //
    // update p, compute stopping criterion
    // dif = ||p-p0||/||p|| = c||b-beta||/||p||
    T nd(0);
    T np(0);
    for (int i =0; i < n; i++) {
      T d = b[i]-beta[i];
      p[i] = p[i] + c*d;
      np = np + p[i]*p[i];;
      nd = nd + d*d;
    }
    np = sqrt(np);
    nd = sqrt(nd);
    dif = fabs(c*nd/(np+1e-8));
    if ( dif < tol )
      break;
  }
  //printf("(%g)\n",dif);
  
  //
  // solution is taken to be the combination of both
  // variables b and beta so that we get more exact zeroes
  //
  for (int i = 0; i < n; i++) {
    if (beta[i] == T(0)) {
      b[i] = T(0);
    }
  }
  return  (j < max_iter);
}
