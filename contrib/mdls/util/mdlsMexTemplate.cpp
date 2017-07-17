#include <mex.h>
#include <math.h>
#ifndef MATLAB_MEX_FILE
  #define MATLAB_MEX_FILE
#endif

/*
 * Input:
 * 0 X ............ input
 *
 * output:
 * 0 Y ............ output
 */

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {

  /* pointers to raw data */
  /* input */
  double *X;
  /* output */
  double *Y;

   /* auxiliary */
  register int i,j;
  int m,n;

   if (nrhs < 1) {
     mexErrMsgTxt("Short help goes here.\n");
       return;
   }

   /* input */
   X = mxGetPr(prhs[0]);   
   m = (int)mxGetM(prhs[0]);
   n = (int)mxGetN(prhs[0]);

   /* initialize output */
   plhs[0]=mxCreateDoubleMatrix(m,n,mxREAL);
   Y=mxGetPr(plhs[0]);

   /* sample loop over all elements */
   /* faster if done 1-D (for element wise operations) */
   /* matlab gives raw data in column-major order, that is
    * columns are concatenated into one single linear vector */
   for ( i = 0 ; i < m ; i++ ) {
     for ( j = 0 ; j < n ; j++ ) {
       /* do something here */
     }
   } 
}
