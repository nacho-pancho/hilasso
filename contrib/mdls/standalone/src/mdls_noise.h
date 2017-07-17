#ifndef MDLS_NOISE
#define MDLS_NOISE

#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include "rand48.h"

#define MRAND48_MAX 0x7FFFFFFFL

template<class T> 
void add_impulse_noise(T* x, T* z, size_t nrows, size_t ncols, double delta, int amax)
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

template<class T> 
void add_symmetric_noise(T* x, T* z, size_t nrows, size_t ncols, double delta, int amax)
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
#endif
