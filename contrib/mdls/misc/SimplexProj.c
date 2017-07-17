/*
 * Matlab C MEX function for projecting a vector to a Simplex
 *
 * Author  : Ignacio Ramirez  <nacho@fing.edu.uy>
 *
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mex.h"

#define DEBUG 0
#define MRAND48_MAX 0x7FFFFFFFL

int comp_double(const void* pa, const void* pb) {
  double a,b;
  
  a = *((double *)pa);
  b = *((double *)pb);
  
  if (a < b) {
    return -1;
  } else if (a > b) {
    return 1;
  } else {
    return 0;
  }
}

void simplex_proj(const double* v, int m, double z, double* w) {
  int t;
  double u;
  memcpy( w, v, m*sizeof(double) );
  /* 
   * 1. sort w
   */
  qsort( w, m, sizeof(double), comp_double );
  /*
   * 2. find theta
   */
  t = 1;
  u = w[m-t]-z;
  while ( ( ((double)t)*w[m-t] ) > u ) {
    t++;
    u += w[m-t];
  }
  u = u / ( (double)t );
  /* 
   * 3. compute solution
   */
  for (t = 0; t < m; t++) {
    w[t] = (v[t] > u) ? v[t]-u : 0.0;
  }
}

/**
 * Z=SimplexProj(X,type,params)
 *
 * inputs:
 *
 * 0 v ....... vector to project. If this is a matrix the procedure will be carried on each column
 * 1 z ....... simplex to project into is w>=0, sum(w)=z
 * 
 * w ......... projection, same dimension as v
 */
void mexFunction(int nlhs,mxArray *plhs[],
		 int nrhs,const mxArray *prhs[])

{
  double *v,*w;
  double z;  
  int n,m;
  int i,j;

  if (nrhs < 1){
    mexErrMsgTxt("Usage: w=SimplexProj(v,z), z optional, defaults to 1.");
    return;
  }
  /*
   * initialize input and output structures
   */
  v = (double*) mxGetPr(prhs[0]);
  n = mxGetN(prhs[0]); 
  m = mxGetM(prhs[0]);
  z = 1.0;
  if (nrhs > 1) {
    z = *mxGetPr(prhs[1]);
  }

  plhs[0] = mxCreateNumericMatrix(m,n,mxDOUBLE_CLASS, mxREAL);
  w = mxGetPr(plhs[0]);

  if ((j == 1) || (i == 1)) {
    simplex_proj(v, m, z, w);
  }

  for (j=0; j < n; j++) {
    simplex_proj(v+j*m, m, z, w+j*m);
  }  
}	
