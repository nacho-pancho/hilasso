#include <mex.h>
#include <math.h>
#include <iostream>

#ifndef MATLAB_MEX_FILE
  #define MATLAB_MEX_FILE
#endif

/*
 * 
 * Perform the canonical transformation of a patch or set of patches.
 * 
 * The transformation is defined so that the brightest corner of the patch
 * is located at the upper-left corner, and the next brightest one is in the
 * upper-right.
 *
 * Input:
 * 0 X ............ original patches
 * 1 TM ........... transformations array (0-offset)
 * 2 IT ........... INVERSE transform: if specified, this vector indicates the transforms that
 *                  were used for each patch, and the inverse are applied to reconver the original
 *                  patches.
 *
 * output:
 * 0 Y ............ transformed patches
 * 1 T ............ index to transformation applied to each patch (0-offset)
 */

typedef enum {
  TRANS_00R=0,
  TRANS_90R=1,
  TRANS_90=2,
  TRANS_270R=3,
  TRANS_270=4,
  TRANS_180R=5,
  TRANS_180=6,
  TRANS_00=7} trans_t;

/* function to determine which transformation corresponds to the patch */
template<class T>  trans_t find_transform(const T* x, const unsigned int w);

/* function to determine which transformation corresponds to the patch */
template<class T> void do_transform(const T* x, const unsigned int m, const short *TM, const trans_t t, T* y);

/* function to determine which transformation corresponds to the patch */
template<class T> void do_transform_inv(const T* x, const unsigned int m, const short *TM, const trans_t t, T* y);

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {

  /* pointers to raw data */
  /* input */
  double *X;
  short *TM;
  short *IT;

  /* output */
  double *Y;
  double *x,*y;
  short *TI;

   /* auxiliary */
  register unsigned int j;
  unsigned int m,n;
  unsigned short w;
  trans_t t;
  int d;

   if (nrhs < 2) {
     mexErrMsgTxt("Short help goes here.\n");
       return;
   }

   /* input */
   X = mxGetPr(prhs[0]);   
   m = (unsigned int)mxGetM(prhs[0]);
   n = (unsigned int)mxGetN(prhs[0]);
   w = (unsigned short)sqrt(m);
   
   TM = (short*) mxGetPr(prhs[1]);
   if ((mxGetM(prhs[1]) !=  w*w) || (mxGetN(prhs[1]) != 7)) {
     mexErrMsgTxt("Invalid transform matrix.");
     return;
   }
   d = 0;
   IT = NULL;
   if (nrhs == 3) {
     IT = (short*) mxGetPr(prhs[2]);
     if (mxGetN(prhs[2]) != n) {
       mexErrMsgTxt("Inverse mode: invalid transform indexes.");
       return;
     }
   }

   /* initialize output */
   plhs[0]=mxCreateDoubleMatrix(m,n,mxREAL);
   Y=mxGetPr(plhs[0]);
   plhs[1]=mxCreateNumericMatrix(1,n,mxINT16_CLASS,mxREAL);
   TI=(short*) mxGetPr(plhs[1]);
   if (IT==NULL) { // direct
     for ( j = 0 ; j < n ; j++ ) {     
       x = &X[m*j];
       y = &Y[m*j];
       t = find_transform(x,w);
       do_transform(x,m,TM,t,y);
       TI[j]=t;
     }
   } else { //inv
     for ( j = 0 ; j < n ; j++ ) {     
       x = &X[m*j];
       y = &Y[m*j];
       t = (trans_t) IT[j];
       do_transform_inv(x,m,TM,t,y);
       TI[j]=t;
     }
   }
}

enum {QUAD_UL=0,QUAD_UR=1,QUAD_DL=2,QUAD_DR=3};

template<class T>  trans_t find_transform(const T* patch, const unsigned int w)
{
  //
  // 1. compute sums of the four  quadrants
  //
  T sum_dl,sum_ul,sum_ur,sum_dr;
  T sum_max;
  int quad_max;
  register unsigned int i,j;
  unsigned int w2;
  trans_t t;
  
  w2 = w/2;
  sum_ul=0;
  sum_dl=0;
  sum_ur=0;
  sum_dr=0;
  for (j=0; j < w2; j++) {
    for (i=0; i < w2; i++) {
      sum_ul += patch[j*w+i];
    }    
    for (i=w-w2; i < w; i++) {
      sum_dl += patch[j*w+i];      
    }    
  }
  for (j=w-w2; j < w; j++) {
    for (i=0; i < w2; i++) {
      sum_ur += patch[j*w+i];
    }    
    for (i=w-w2; i < w; i++) {
      sum_dr += patch[j*w+i];      
    }    
  }
#ifdef DEBUG
  std::cerr << "UL=" << sum_ul << " UR=" << sum_ur <<  " DL=" << sum_dl << " DR=" << sum_dr << std::endl;
#endif
  //
  // 2. determine transformation
  // 
  quad_max = QUAD_UL;
  sum_max = -1e20;
  if (sum_ur > sum_max) {
      quad_max = QUAD_UR;
      sum_max = sum_ur;
  }
  if (sum_dl > sum_max) {
      quad_max = QUAD_DL;
      sum_max = sum_dl;
  }
  if (sum_dr > sum_max) {
      quad_max = QUAD_DR;
      sum_max = sum_dr;
  }
  switch (quad_max)
    {
    case QUAD_UL:
      t = (sum_dl > sum_ur)?  TRANS_00R: TRANS_00;
     break;
    case QUAD_UR:
      t = (sum_ul > sum_dr)? TRANS_90R: TRANS_90;
      break;
    case QUAD_DL:
      t = (sum_dr > sum_ul)? TRANS_270R : TRANS_270;
      break;
    case QUAD_DR:
      t = (sum_ur > sum_dl)? TRANS_180R : TRANS_180;
      break;
    }  
#ifdef DEBUG
  std::cerr << "quad_max=" << quad_max << " sum_max=" << sum_max <<  " t=" << (int)t << std::endl;
#endif
  return t;
}

template<class T> void do_transform(const T* x, const unsigned int m, const short *transforms, const trans_t t, T* y)
{
  register unsigned int i;
  const short* trans;
#ifdef DEBUG
  std::cerr << "m=" << m << " t=" << t << std::endl;
#endif
  if (t == TRANS_00) {
    for (i=0; i < m; i++) {
      y[i]=x[i];
    }    
  } else {
    trans = &transforms[m*t];
    for (i=0; i < m; i++) {
      y[trans[i]] = x[i];
    }
  }
}

template<class T> void do_transform_inv(const T* x, const unsigned int m, const short *transforms, const trans_t t, T* y)
{
  register unsigned int i;
  const short* trans;
#ifdef DEBUG
  std::cerr << "m=" << m << " t=" << t << std::endl;
#endif
  if (t == TRANS_00) {
    for (i=0; i < m; i++) {
      y[i]=x[i];
    }    
  } else {
    trans = &transforms[m*t];
    for (i=0; i < m; i++) {
      y[i] = x[trans[i]];
    }
  }
}
