#include <mex.h>
#include <math.h>
#ifndef MATLAB_MEX_FILE
  #define MATLAB_MEX_FILE
#endif

#include "mex_traits.h"

template<class T> void actual_mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]);

/*
 * Input:
 * 0 I ............ image to be decomposed
 * 1 w ............ width of patches
 * 2 ov ........... overlap (both vert. and horiz)
 *
 * output:
 * 0 X ............ Patches of the image given as column vectors
 * 0 g ............ (2xn) grid coordinates: each column contains the 
 *                  coordinates of a point in the grid from where the
 *                  pixels where taken.
 */
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
  mxClassID id;
  if (nrhs < 1) {
    mexErrMsgTxt("Need at least one argument: image.\n");
  }
  id = mxGetClassID(prhs[0]);
  switch (id) {
    //  case mxCHAR_CLASS:
    //actual_mexFunction<char>(nlhs,plhs,nrhs,prhs);
    //return;
  case mxUINT8_CLASS:
    actual_mexFunction<unsigned char>(nlhs,plhs,nrhs,prhs);
    return;
    //case mxINT16_CLASS:
    //actual_mexFunction<short>(nlhs,plhs,nrhs,prhs);
    //return;
  case mxUINT16_CLASS:
    actual_mexFunction<unsigned short>(nlhs,plhs,nrhs,prhs);
    return;
  case mxINT32_CLASS:
    actual_mexFunction<int>(nlhs,plhs,nrhs,prhs);
    return;
  case mxUINT32_CLASS:
    actual_mexFunction<unsigned int>(nlhs,plhs,nrhs,prhs);
    return;
    //case mxINT64_CLASS:
    //actual_mexFunction<long>(nlhs,plhs,nrhs,prhs);
    //return;
    //case mxUINT64_CLASS:
    // actual_mexFunction<unsigned long>(nlhs,plhs,nrhs,prhs);
    //return;
  case mxSINGLE_CLASS:
    actual_mexFunction<float>(nlhs,plhs,nrhs,prhs);
    return;
  case mxDOUBLE_CLASS:
      actual_mexFunction<double>(nlhs,plhs,nrhs,prhs);
      return;
  default:
    mexErrMsgTxt("Input precision not supported.");
  }
}

template<class T>
void actual_mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
  /* input data */
  T *I;
  int m,n,m2,n2;
  int w,ov;
  int step;

  /* output data */
  T *X;
  int D,N;
  int d;


   /* auxiliary */
  register int i,j,i2,j2;

   if (nrhs < 1) {
     fprintf(stderr,"Need at least one argument: image.\n");
       return;
   }

   /* input data */
   I = (T*) mxGetPr(prhs[0]);
   m = (int)mxGetM(prhs[0]);
   n = (int)mxGetN(prhs[0]);
   if (nrhs > 1)
     w = (int) *mxGetPr(prhs[1]);
   else
     w = 8;

   if (nrhs > 2)
     ov = (int) *mxGetPr(prhs[2]);
   else
     ov = w-1;
   
   D = w*w; /* dimension of patches */
   step = w - ov;
   

   /* these are temporary to compute the number of patches */
   j = (n-w)/step + 1; /* case when image breaks exactly */
   i = (m-w)/step + 1;
   /* 
    * if image doesn't break exactly into patches of this size
    * we pad it virtually.
    */
   if (((n-w) % step) != 0) {
     n2 = ((n-w)/step + 1)*step + w;
     j++;
   } else {
     n2 = n;
   }
   if (((m-w) % step) != 0) {
     m2 = ((m-w)/step + 1)*step + w;
     i++;
   } else {
     m2 = m;
   }
   N = i*j; 
   /* mexPrintf("m=%d n=%d m2=%d n2=%d D=%d N=%d w=%d ov=%d step=%d\n", m, n, m2, n2, D, N, w, ov, step); */
   plhs[0]=mxCreateNumericMatrix(D,N,mex_traits<T>::mx_class(),mxREAL);
   X= (T*) mxGetPr(plhs[0]);
   

   d = 0; /* since X is col-major, X 'wraps' to the next column automatically 
	     when d surpasses D */
   for ( i = 0 ; i <= (m2-w) ; i += step ) {
     for ( j = 0 ; j <= (n2-w) ; j += step ) {
       /* add a patch */
       for ( j2 = 0; j2 < w ; j2++ ) {
          for ( i2 = 0 ; i2 < w ; i2++ ) {
	    if ((i+i2) >= m) 
	      X[d++] = T();
	    else if ((j+j2) >= n)
	      X[d++] = T();
	    else {
	      X[d++] = I[ m*(j+j2)+i+i2 ];
	    }
	 }
       } 
     }
   } 
}
