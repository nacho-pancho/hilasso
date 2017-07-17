/*
 * Matlab C MEX implementation of bwlabel
 * recursive algorithm
 *
 * File    : bwlabel5.c
 * Author  : Ignacio Ramirez  <nacho@fing.edu.uy>
 *
 */
#include <math.h>
/*#include <malloc.h>*/
#include <stdio.h>
#include "mex.h"
//#include "matrix.h"

#define LI(ii,jj) (label_image[(jj)*m+(ii)])
#define XI(ii,jj) (bw_image[(jj)*m+(ii)])

#define DEBUG 0

typedef unsigned int label_t;

int n,m;
int mn;
int stack_top=0;
int *stack;

int debug;

double *bw_image;
double *label_image;

void label(int i, int j, double value) {
  /* label current pixel */
  stack_top=0;
  /* printf("starting at (%d,%d)\n",i,j); */
  do {
    LI(i,j)=value;
    if ( (i>0) &&(XI(i-1,j)>0) && (LI(i-1,j)==0) ) {
      /* process west child */
      stack[stack_top++]=i;
      stack[stack_top++]=j;
      i--; 
    } else if ( (j>0) &&  (XI(i,j-1)>0) && (LI(i,j-1)==0) ) {
      /* process north child */
      stack[stack_top++]=i;
      stack[stack_top++]=j;
      j--; 
    } else if ( (i<(m-1)) &&  (XI(i+1,j)>0) && (LI(i+1,j)==0) ) {
      /* process south child */
      stack[stack_top++]=i;
      stack[stack_top++]=j;
      i++; 
    } else if ( (j<(n-1)) && (XI(i,j+1)>0) && (LI(i,j+1)==0) ) {
      /* process east child */
      stack[stack_top++]=i;
      stack[stack_top++]=j;
      j++; 
    } else if (stack_top > 0) {
      j=stack[--stack_top];
      i=stack[--stack_top];
    }
  } while (stack_top);
}

/**
 * area_hist=areahist(label_image)
 */
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
  int i,j;
  int Lmax;

  /* temporary data */


  //mexErrMsgTxt("Does not work!!");
  /* 
   * check for correct # of input variables 
   */
  if (nrhs<1) {
    mexErrMsgTxt("There are exactly 1 input parameters allowed!");
  }
  if (nrhs>1){
    //mexErrMsgTxt("There are exactly 1 input parameters allowed!");
    printf("debug mode\n");
    debug = 1;
  } else {
    debug = 0;
  }
  /* 
   * get parameters
   */
  bw_image = mxGetPr(prhs[0]);
  n  = mxGetN(prhs[0]); 
  m  = mxGetM(prhs[0]); 
  mn = n*m;
  stack=mxCalloc(mn,sizeof(int)); 
  
  /* 
   * create output histogram and initialize it
   */
  plhs[0]     = mxCreateDoubleMatrix(m,n,mxREAL);
  label_image = mxGetPr(plhs[0]);

  Lmax=0.0;
  for (j=0; j < n; j++) {    
    for (i=0; i < m; i++) {
      if ((XI(i,j) > 0.0) &&(LI(i,j)==0)) {
	label(i,j,(double)++Lmax);
      } 
    } /* for each column element */
  } /* for each column */

  mxFree(stack);
}
