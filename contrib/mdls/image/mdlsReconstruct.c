#include <mex.h>
#include <math.h>
#ifndef MATLAB_MEX_FILE
  #define MATLAB_MEX_FILE
#endif

/*
 * Input:
 * 0 X ............ patches (D x N)
 * 1 m ............ height of output image
 * 2 n ............ width of output image
 * 3 overlap ...... overlap of patches
 % 4 DC ........... optional DC, (1 x N)
 *
 * output:
 * 0 I ............ Output image.
 */

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {

   /* data */
  double *X, *I;
  int ov;
  int D,N;
  int m,n,L,LX;
  int d;
  int w;
  int step;
  double *DC = NULL;
  double dc;

   /* auxiliary */
  register int i,j,i2,j2,idx;
   mxArray *mxR;
   double *R;

   if (nrhs < 4) {
     fprintf(stderr,"Need at least four arguments: X, m, n, overlap.\n");
       return;
   }

   /* input data */
   if (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS) {
     mexErrMsgTxt("Currently only double precision X is supported. Sorry!");
   }
   X = mxGetPr(prhs[0]);
   
   /* dimensions */
   D  = (int)mxGetM(prhs[0]);
   N  = (int)mxGetN(prhs[0]);

   w = (int) sqrt(D);

   m = *mxGetPr(prhs[1]);
   n = *mxGetPr(prhs[2]);
   L = m*n;

   ov = (int) *mxGetPr(prhs[3]);
   if (ov >= w) {
     mexErrMsgTxt("Invalid overlap.");
   }
   step = w - ov;
   /*
    * check whether the size of X is consistent with the requested
    * size of the resulting image.
    */
   LX = ((int)ceil( ((double)(m-w)) / ((double)step) ) + 1 ) 
     * ((int)ceil( ((double)(n-w)) / ((double)step) ) + 1);
   if (LX != N)
     mexErrMsgTxt("Size of X is incompatible with m,n,overlap.");
   if (nrhs > 4) {
     DC = mxGetPr(prhs[4]);
   }
   /*
    *  fprintf(stderr,"m=%d n=%d L=%d D=%d N=%d w=%d ov=%d step=%d\n",
    *	   m,n,L,D,N,w,ov,step);
    */
   /* output */
   plhs[0]=mxCreateDoubleMatrix(m,n,mxREAL);
   I=mxGetPr(plhs[0]);

   mxR=mxCreateDoubleMatrix(m,n,mxREAL);
   R=mxGetPr(mxR);

   /* initialize auxiliary */
   for (i = 0; i < L ; i++) {
     R[i] = 0.0; I[i] = 0.0;
   }
   d = 0;
   for ( i = 0 ; i <= (m-w) ; i += step ) {
     for ( j = 0 ; j <= (n-w) ; j += step ) {
       /* d = 0; */
       /* add a patch */
       dc = (DC!=NULL) ? DC[d/D] : 0.0;
       for ( j2 = 0; j2 < w ; j2++ ) {
          for ( i2 = 0 ; i2 < w ; i2++ ) {
	   idx = m*(j+j2)+i+i2;
	   I[idx] += (X[d++] + dc); 
	   R[idx] += 1.0 ;
	 }
       } 
     }
   } 

   /* normalize */
   for (i = 0; i < L ; i++) {
     if (R[i] > 0.0)
       I[i] /= R[i];
   }
   
  /* destroy temporary space */
  mxDestroyArray(mxR);
}
