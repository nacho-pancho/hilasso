#include <mex.h>
#include <math.h>
#ifndef MATLAB_MEX_FILE
  #define MATLAB_MEX_FILE
#endif

/*
 * 
 * Initialize canonical transformation matrix given a patch width.
 *
 * Input:
 * 0 w ............ patch width
 *
 * output:
 * 0 T ............ canonical trans. matrix
 */

#define NUM_TRANSFORMS  7 /* 4 rotations * 2 reflections - identity */

typedef enum {
  TRANS_00R=0,
  TRANS_90R=1,
  TRANS_90=2,
  TRANS_270R=3,
  TRANS_270=4,
  TRANS_180R=5,
  TRANS_180=6,
  TRANS_00=7} trans_t;

/* function that precomputes the transforms for a given patch width */
void initialize_transforms(const int w, short int *T);

/* 
 * each transform is a permutation, so it is stored as an array of indexes 
 * so that Y[i]=X[trans[i]]
 */
short int *transforms[NUM_TRANSFORMS];

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {

  /* pointers to raw data */
  /* input */
  int w;

  /* output */
  short int *T;

   if (nrhs < 1) {
     fprintf(stderr,"Short help goes here.\n");
       return;
   }

   /* input */
   w = (int) *mxGetPr(prhs[0]);   

   /* 
    * Initialize output. The transf. matrix has 7 columns, each one corresponding 
    * to a permutation of the indexes of original to permuted patches 
    */
   plhs[0]=mxCreateNumericMatrix(w*w,7,mxINT16_CLASS,mxREAL);
   T=(short int*) mxGetPr(plhs[0]);

   /* sample loop over all elements */
   /* faster if done 1-D (for element wise operations) */
   /* matlab gives raw data in column-major order, that is
    * columns are concatenated into one single linear vector */
   initialize_transforms(w,T);
}


void initialize_transforms(const int w, short int *T)
{
  int i,j,jt,it;
  int n;
  n = w*w;
  int r;
  for (i=0; i < w; i++) {
    for (j=0; j < w; j++) {
      r = j*w+i; /* matlab indexes are from 1 */
      // 0: 0 rotaation + flip      
      it = j; 
      jt = i; 
      T[    jt*w+it] = r;
      // 1: +90 rotation  + flip along diagonal
      it = i; 
      jt = w-1-j; 
      T[  n+jt*w+it] = r;
      // 2: +90 rotation only
      it = w-1-j; 
      jt = i; 
      T[2*n+jt*w+it] = r;
      // 3: -90 rotation + flip
      it = w-1-i;
      jt = j; // rot+flip
      T[3*n+jt*w+it] = r;
      // 4: -90 rotation only
      it = j;
      jt = w-1-i; // rot only
      T[4*n+jt*w+it] = r;
      // 5: 180 rotation + flip
      it = w-1-j;
      jt = w-1-i; //rot+flip
      T[5*n+jt*w+it] = r;
      // 6: 180 rotation only
      it = w-1-i;
      jt = w-1-j; // rot only
      T[6*n+jt*w+it] = r;
    }
  }
}
