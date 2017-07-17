/*
 * Matlab C MEX function for adding noise to some data
 *
 * Author  : Ignacio Ramirez  <nacho@fing.edu.uy>
 *
*/

#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include "mex.h"
#include "rand48.h"

#define MRAND48_MAX 0x7FFFFFFFL

template<class T> void add_impulse_noise(T* x, T* z, size_t nrows, size_t ncols, double delta, int amax)
{
  size_t n;
  n=nrows*ncols;
  double ps,pp,p;
  pp=delta/2.0;
  ps=1.0-pp;
  for(size_t i=0; i < n; i++) {
    p = drand48();
    if ( p < pp)
      z[i]=0;
    else if (p < ps)
      z[i]=x[i];
    else
      z[i]=amax;
  }
}

template<class T> void add_symmetric_noise(T* x, T* z, size_t nrows, size_t ncols, double delta, int amax)
{
  size_t n;
  n=nrows*ncols;
  double p;
  double pes;
  pes = ((double)amax-1.0)/delta;
  for(size_t i=0; i < n; i++) {
    p = drand48();
    if ( p < delta ) {  /* error */
      /*
       * divide range [0,perr) evenly among the rest of the asize-1 symbols
       */
      z[i] = p*pes;
      if (z[i] >= x[i]) /* skip original value, as this is considered in the 'else' clause */
        z[i]++;
    }
    else { /* clean */
      z[i] = x[i];
    }
  }
}

template<class T> void add_gaussian_noise(T* x, T* z, size_t nrows, size_t ncols, double sigma, int amax)
{
  size_t n;
  double v, w, p, q, r2,y;
  double tol;

  n=nrows*ncols;
  tol = 1e-10; // added by nacho to avoid numerical problems
  for(size_t i=0; i < n; i++) {
    /* p = drand48(); */
    /* Box-Muller transform:
     * Choose v,w in uniform square (-1,-1) to (+1,+1)
     * until (v,w) is in the unit circle and different from (0,0)
     */
    do
    {
      p= ((double) lrand48()) / ( (double)MRAND48_MAX + 1.0 );
      q= ((double) lrand48()) / ( (double)MRAND48_MAX + 1.0 );
      v = -1.0 + 2.0 * p;
      w = -1.0 + 2.0 * q;
      r2 = v * v + w * w;
    } while ( ((1.0-r2) < tol) || ((r2-0.0) < tol) ); 
    y= ((double)x[i]) + sigma * w * sqrt (-2.0 * log (r2) / r2);
    if (y > (double)amax)
      z[i]=amax;
    else if (y < 0.0)
      z[i]=0;
    else
      z[i]=(int)(y+0.5);
  }
}

/**
 * Z=addnoise(X,type,params)
 *
 * inputs:
 *
 * 0 Z ....... clean image
 * 1 type .... type of noise. Can be 'Gaussian','Impulse','Symmetric'
 * 2 params .. parameters of noise. Depends on type of noise.
 * 
 * Z ......... noisy image
 * stats ..... if specified, return the histogram
 *
 */
void mexFunction(int nlhs,mxArray *plhs[],
		 int nrhs,const mxArray *prhs[])

{
  unsigned char *x,*z;
  double *params;  
  int n,m;
  int seed;

  mxChar* type_mx;
  char type[16];

  if (nrhs < 3){
    mexErrMsgTxt("At least three parameters are required: clean image, noise type, and parameters. A fourth parameter is the seed to be used for the RNG.");
    return;
  }
  if (nrhs >= 4){
    seed = (int) *mxGetPr(prhs[3]);
  } else {
    seed = (int) time(NULL);
  }
  srand48(seed);

  if (!mxIsUint8(prhs[0])) {
    mexErrMsgTxt("Input data must be uint8");
  }
  /*
   * initialize input and output structures
   */
  x = (unsigned char*) mxGetPr(prhs[0]);
  n = mxGetN(prhs[0]); 
  m = mxGetM(prhs[0]);

  if (mxGetChars(prhs[1]) == NULL)
  { /* we didn't receive a char string here */
    mexErrMsgTxt("Second parameter must be a string");
  }
  mxGetString(prhs[1],&type[0],16);
  params = mxGetPr(prhs[2]);
  plhs[0] = mxCreateNumericMatrix(m,n,mxUINT8_CLASS, mxREAL);
  z = (unsigned char*) mxGetPr(plhs[0]);

  if (!strcmp(type,"Gaussian"))
  {
    add_gaussian_noise<unsigned char>(x,z,m,n,*params,255);
  } else if (!strcmp(type,"Impulse"))
  {
    add_impulse_noise<unsigned char>(x,z,m,n,*params,255);
  } else if (!strcmp(type,"Symmetric"))
  {
    add_symmetric_noise<unsigned char>(x,z,m,n,*params,255);
  } 
}	
