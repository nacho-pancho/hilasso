/*
 * Matlab C MEX function for adding noise to some data
 *
 * Author  : Ignacio Ramirez  <nacho@fing.edu.uy>
 *
*/

#include <math.h>
#include <stdio.h>
#include "mex.h"

#define DEBUG 0

/**
 * Y=SoftThresholding(W,T)
 *
 * inputs:
 *
 * 0 X ....... input values
 * 1 T ....... threshold, may be same size as X or scalar
 *
 * outputs:
 * Y ......... thresholded data
 * d ......... cumulative absolute difference between input and output (optional)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *X, *T, *Y;
  double t;
  double d;

  int n,m,N,i;
  double nt,mt;
  
  char* type;
  if (nrhs != 2){
    mexErrMsgTxt("Two parameters are required. See help.");
    return;
  }
  /*
   * initialize input and output structures
   */
  X = mxGetPr(prhs[0]);
  n = mxGetN(prhs[0]); 
  m = mxGetM(prhs[0]);
  N=n*m;

  T = mxGetPr(prhs[1]);
  nt = mxGetN(prhs[1]); 
  mt = mxGetM(prhs[1]);

  if (nt*mt > 1) {  
    if ((nt != n) ||  (mt != m)) {
	    mexErrMsgTxt("Threshold must be same size as X or scalar.");
    	return;
    }
	} 
  plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);

  Y = mxGetPr(plhs[0]);
  d = 0;
  
  if (nt*mt > 1) {  
    /* T is a matrix */
    for (i=0; i < N; i++) {
      t=T[i];
      if (X[i] < -t) {
        Y[i] = X[i] + t ; d += t;
      } else if (X[i] > t) {
        Y[i] = X[i] - t ; d += t;
      } else {
        Y[i]=0.0; 
      }
    }
  }  else { /* T is a scalar */
    t = T[0];
    for (i=0; i < N; i++) {
      if (X[i] < -t) {
        Y[i] = X[i] + t ; d += t;
      } else if (X[i] > t) {
        Y[i] = X[i] - t ; d += t;
      } else {
        Y[i]=0.0;
      }
    }
  }

  if (nlhs == 2) {
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(plhs[1]) = d;
  }

}
