/*
 * Matlab C MEX implementation of PSNR and MSE measures
 *
 * File    : imgdif.c
 * Author  : Ignacio Ramirez  <nacho@fing.edu.uy>
*/

#include <math.h>
#include <stdio.h>
#include "mex.h"
#include <stdlib.h>

#define DEBUG 0

/**
 * R=imgdif(X,Z)
 *
 * inputs (prhs):
 * 0 X: clean image
 * 1 Z: noisy image
 *
 * outputs (plhs):
 * 0 PSNR
 * 1 MSE
 */
void mexFunction(int nlhs,mxArray *plhs[],
		 int nrhs,const mxArray *prhs[])

{
  /* data */
  double *X,*Z;
  /* sizes and ranges */
  int N,M,L;
  /* indexes */
  register int n;
  /* other */
  register int d;
  double tmp;
  double mse;
  double  *pPSNR, *pMSE;

  /* 
   * check for correct # of input variables 
   */
  if (nrhs!=2)
  {
    mexErrMsgTxt("There are exactly 2 input parameters!");
    return;
  } 
  
  /* 
   * initialize input and output structures
   */
  X = mxGetPr(prhs[0]);
  N = mxGetN(prhs[0]); 
  M = mxGetM(prhs[0]);
  L = M*N;
  Z = mxGetPr(prhs[1]);

  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  if (nlhs > 1)
  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);

  pPSNR = mxGetPr(plhs[0]);
  if (nlhs > 1)
    pMSE = mxGetPr(plhs[1]);
  
  tmp=0.0;
  for (n=0; n<L; n++)
  {
    d = X[n]-Z[n];
    tmp += d*d; 
  }
  mse = tmp/((double)L);
  *pPSNR = 10.0*log10(255.0*255.0/ mse);
  if (nlhs > 1)
    *pMSE = mse;
}
