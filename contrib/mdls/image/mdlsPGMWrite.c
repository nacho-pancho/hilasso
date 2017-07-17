	/*
	 * Read 8 or 16 bit PGM files
	 *
	 * File    : pgmread.c
	 * Author  : Ignacio Ramirez  <nacho@fing.edu.uy>
	 *
	 */

#include <stdio.h>
#include "mex.h"
#include <stdlib.h>

typedef unsigned char buffer_t;

/**
 * I=pgmread(fname)
 *
 * inputs (prhs):
 * 0 fname: image name
 *
 * outputs (plhs):
 * 0 I : image
 */
void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  char fname[1024];
  FILE *f = NULL;
  double *I;
  int M, N, MN;
  size_t N2;
  size_t i,j,k;
  int lsb,msb;
  unsigned char buffer[4096];
  double minv,maxv,scale;

  /*
   * check for correct # of input variables
   */
  if (nrhs != 2)
    {
      mexErrMsgTxt ("There are 2 input parameters allowed!");
      return;
    }
  M = N = 0;
  mxGetString (prhs[0], fname, 1024);
  f = fopen (fname, "w");
  if (f == NULL)
    {
      printf ("Error opening %s for writing\n", fname);
      return;
    }
  I = mxGetPr (prhs[1]);
  M = mxGetM(prhs[1]);
  N = mxGetN(prhs[1]);
  MN=M*N;
  minv=1e15;
  maxv=-1e15;
  for (i=0; i < MN; i++) {    
    if (I[i] > maxv)
      maxv=I[i];
    else if (I[i] < minv)
      minv=I[i];
  }
  if ((maxv < 256.0) && (minv >=0)) {
    scale=1.0;
    minv=0.0;
  } else {
    scale=255.0/(maxv-minv);
  }
  
  fprintf(f,"P5 %d %d 255\n",M,N);
  k=0;
  for (i=0; i < M; i++) {
    for (j=0; j < N; j++) {
      buffer[k++]=(unsigned char)((I[i+j*N]-minv)*scale +0.5);
      if (k >= 4096) {
        fwrite(buffer,1,sizeof(unsigned char)*4096,f);
        k=0;
      }
    }
  }
  if (k > 0) {
    fwrite(buffer,1,sizeof(unsigned char)*k,f);
   }
  fclose(f);
}
