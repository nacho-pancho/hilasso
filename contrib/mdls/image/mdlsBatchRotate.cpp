#include <mex.h>
#include <math.h>
#ifndef MATLAB_MEX_FILE
  #define MATLAB_MEX_FILE
#endif

#include "mex_traits.h"

template<class T> void actual_mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]);

/*
 * Rotation of a sqrt(n)xsqrt(n) patch is represented as a linear mapping from R^n -> R^n. Y = H*X
 * Each column in H represents the mapping of a single pixel from X onto Y.
 *
 * Input:
 * 0 X ............ (mxn) images to be rotated, given as columns of a matrix
 * 1 R ............ (n^2xr) transformation matrices. Each column corresponds to a square transformation matrix.
 * 2 J ............ (1xn) vector of transformation indexes. Each element J(i) is between 1 and r, indicating
 *                  which transformation is to be applied to the corresponding column in X.
 *
 * output:
 * 0 Y ............ Rotated (transformed) images.
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
    //case mxUINT8_CLASS:
    //actual_mexFunction<unsigned char>(nlhs,plhs,nrhs,prhs);
    //return;
    //case mxINT16_CLASS:
    //actual_mexFunction<short>(nlhs,plhs,nrhs,prhs);
    //return;
    //case mxUINT16_CLASS:
    //actual_mexFunction<unsigned short>(nlhs,plhs,nrhs,prhs);
    //return;
    //case mxINT32_CLASS:
    //actual_mexFunction<int>(nlhs,plhs,nrhs,prhs);
    //return;
    //case mxUINT32_CLASS:
    //actual_mexFunction<unsigned int>(nlhs,plhs,nrhs,prhs);
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
  T *X, *R, *J;
  int m, n, m2, r;
  T* Rj, *Xj, *Yj;
  T t;

  /* output data */
  //typedef typename mex_traits<T>::op_t out_t;
  T* Y; 

   /* auxiliary */
  register int i,j,k;

   if (nrhs != 1) {
     mexErrMsgTxt("Need exactly three arguments.\n");
     return;
   }

   /* 
    * INPUT DATA
    */
   X = (T*) mxGetPr(prhs[0]);
   m = (int)mxGetM(prhs[0]);
   n = (int)mxGetN(prhs[0]);

   R = (T*) mxGetPr(prhs[1]);
   m2 = (int)mxGetM(prhs[1]);
   r = (int)mxGetN(prhs[1]);
   if ( ((int)sqrt(m2)) != m) {
     mexErrMsgTxt("Invalid rotation matrix.\n");
     return;
   }
   
   J = (T*) mxGetPr(prhs[2]);
   if ( (mxGetN(prhs[2]) != n)  && (mxGetM(prhs[2]) != m) ) {
     mexErrMsgTxt("Invalid rotation indexes vector.\n");
     return;
   }

   /*
    * INITIALIZE OUTPUT
    */
   plhs[0] = mxCreateNumericMatrix(m,n,mex_traits< T >::mx_class(),mxREAL);
   Y = (T*) mxGetPr(plhs[0]);

   /*
    * PROCESS
    */
   for (j=0; j < n; j++) {
     /* select transformation */
     int ridx = (int) (J[j]-1);
     Rj = &R[m2*ridx];
     Xj = &X[m*j];
     Yj = &Y[m*j]; /* Matlab already initializes this to 0 */
     /* matrix-vector mult., col major: this is, 
      * add each column of Rj one at a time, multiplied by Xj[k] */
     for (k=0; k < m; k++) {       
       t = Xj[k];
       for (i=0; i < m; i++) {
	 Yj[i] += t * *(Rj++);
       }
     }
   }
}
