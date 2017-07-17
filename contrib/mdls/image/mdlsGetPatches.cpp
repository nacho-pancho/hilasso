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
 * 2 grid ......... (2xN) [x;y] center of the patches in image (indexes start at 0)
 *                  for patches with even width, center is (w/2+1,w/2+1) pixel in patch
 * output:
 * 0 X ............ Patches of the image given as column vectors
 * 1 DC ........... DC of patches (optional)
 * 2 VAR .......... Variance of patches (optional)
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

#define MACRO_clip_coords(x,n) (  ( (((x) >=0) ? (x) : 0 ) < (n) )? (x) : ((n)-1)  )

template<class T>
void actual_mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
  /* input data */
  T *I;
  int m,n,m2,n2;
  int w,w2;
  int step;
  int* grid;
  /* output data */
  typedef typename mex_traits<T>::op_t out_t;
  out_t* X; 
  double *DC = NULL,*VAR = NULL;
  int D,N;
  double dc,var;
  int d;
  bool computeDC,computeVAR;
  
   /* auxiliary */
  register int i,j,k,i2,j2,i3,j3;

   if (nrhs < 3) {
     mexErrMsgTxt("Needs three arguments: image, width and grid.\n");
       return;
   }

   /* input data */
   I = (T*) mxGetPr(prhs[0]);
   m = (int)mxGetM(prhs[0]);
   n = (int)mxGetN(prhs[0]);
   w = (int) *mxGetPr(prhs[1]);
   D = w*w; /* dimension of patches */
   w2  = w/2; /* patch radius rounded to nearest smallest integer */

   if (mxGetClassID(prhs[2]) != mxINT32_CLASS) {
     mexErrMsgTxt("grid needs to be INT32. Hint: use INT32(grid).");
     return;
   }

   grid = (int*) mxGetPr(prhs[2]);
   if (mxGetM(prhs[2]) != 2) {
     mexErrMsgTxt("grid needs to be a 2xN matrix where each column is an [x;y] pair.");
     return;
   }
   N = mxGetN(prhs[2]);

   plhs[0] = mxCreateNumericMatrix(D,N,mex_traits< out_t >::mx_class(),mxREAL);
   X = (out_t*) mxGetPr(plhs[0]);

   if (nlhs > 1) {
     plhs[1] = mxCreateNumericMatrix(1,N,mxDOUBLE_CLASS,mxREAL);
     DC = (double*)mxGetPr(plhs[1]);
     computeDC = true;
   } else {
     computeDC = false;
   }

   if (nlhs > 2) {
     plhs[2] = mxCreateNumericMatrix(1,N,mxDOUBLE_CLASS,mxREAL);
     VAR = (double*)mxGetPr(plhs[2]);
     computeVAR = true;
   } else {
     computeVAR = false;
   }

   d = 0; /* since X is col-major, X 'wraps' to the next column automatically 
	     when d surpasses D */
   for (k = 0; k < (2*N); ) {
     j = grid[k++]-w2; // j indexes the 'x' direction
     i = grid[k++]-w2; // i indexes the 'y' dir
     /* add a patch */
     dc = 0.0;
     var = 0.0;
     for ( j2 = 0; j2 < w ; j2++ ) {
       for ( i2 = 0 ; i2 < w ; i2++ ) {
	 T x;
	 i3 = MACRO_clip_coords(i+i2,m);
	 j3 = MACRO_clip_coords(j+j2,n);
	 X[d++] = x = I[ m*j3 + i3 ];	 
	 dc += x;
	 var += x*x;
       } // i2
     } // j2, inner loop to copy patch data
     if (computeDC) {
       /* compute DC */
       dc /= double(D);
       DC[d/D-1] = dc;
       /* remove DC from patch. d is now on the first element of the next patch */
       //for (j2 = 1; j2 <= D; j2++)
       // X[d-j2] -= T(dc);
       if (computeVAR) {
	 var = var/double(D) - dc*dc; // BIASED estimator
	 VAR[d/D-1] = var;
       }
     }
   } // for each  grid location
}
