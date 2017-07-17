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

int read_pnm_header (FILE * infile, int *colsp, int *rowsp, int *maxsp,
		     int *typep);
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
  buffer_t* buf=NULL;
  FILE *f = NULL;
  double *I;
  int M, N, maxs, type;
  size_t N2;
  size_t i,j;
  int lsb,msb;
  
  /*
   * check for correct # of input variables
   */
  if (nrhs != 1)
    {
      mexErrMsgTxt ("There is 1 input parameters allowed!");
      return;
    }
  M = N = 0;
  mxGetString (prhs[0], fname, 1024);
  f = fopen (fname, "rb");
  if (f == NULL)
    {
      printf ("Error opening %s for reading\n", fname);
      plhs[0] = mxCreateDoubleMatrix (M, N, mxREAL);
      return;
    }
  read_pnm_header(f, &N,&M,&maxs,&type);
  plhs[0] = mxCreateDoubleMatrix (M, N, mxREAL);
  I = mxGetPr (plhs[0]);
  if (maxs < 256) { /* 16 bit values */
	buf = (buffer_t*) malloc(sizeof(buffer_t)*N);
	for (i = 0; i<M; i++)
	{
		fread(buf,sizeof(unsigned char),N,f);
		for (j = 0; j < N; j++)
  			I[M*j+i] = buf[j];
	}
  } 
  else 
    { /* 16 bit values */
    buf = (buffer_t*) malloc(sizeof(buffer_t)*N*2);
    for (i = 0; i<M; i++)
    { 
      fread(buf,sizeof(char),2*N,f);
      for (j = 0; j < N; j++) {
        /* MSB then LSB */
	msb = buf[2*j];
	lsb = buf[2*j+1];
        I[M*j+i] = ((msb<<8)&0xff00) | (lsb&0xff);
      }
     }
  }
  free(buf);
  fclose(f);
}

#define MAXLINE 1024
#define LINELEN MAXLINE

int
read_pnm_header (FILE * infile, int *colsp, int *rowsp, int *maxsp,
		 int *typep)
{
  int rows=-1, cols=-1, maxs=-1, type=-1;
  int res;
  char line[LINELEN + 1];
 
  while (fgets (line, LINELEN, infile) != NULL
	 && (line[0] == '#' || line[0] == 0 || line[0] == '\n'));

  res = sscanf (line, "P%d %d %d %d", &type,&cols,&rows,&maxs);
  if ((res < 1) || (type != 5 && type != 6 && type != 7))
    {
      fprintf (stderr, "bad header: type\n");
      return -1;
    }
  if (res < 4) { /* read header, need cols rows and maxs */
    while (fgets (line, LINELEN, infile) != NULL
  	 && (line[0] == '#' || line[0] == 0 || line[0] == '\n'));
    res = sscanf (line, "%d %d %d", &cols, &rows,&maxs);
    if (res < 2)
    {
      fprintf (stderr, "bad header: cols rows\n");
      return -1;
    }
    if (res < 3)  { /* read header, cols, rows; need maxs */
      while (fgets (line, LINELEN, infile) != NULL
  	 && (line[0] == '#' || line[0] == 0 || line[0] == '\n'));
      if (sscanf (line, "%d", &maxs) != 1)
      {
        fprintf (stderr, "bad header: maxs\n");
        return -1;
      }
    }
  } 
  *colsp = cols;
  *rowsp = rows;
  *maxsp = maxs;
  *typep = type;
  return 0;
}

