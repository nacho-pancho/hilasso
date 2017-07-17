/*
 * Matlab C MEX implementation of area histogram
 *
 * File    : areahist.c
 * Author  : Ignacio Ramirez  <nacho@fing.edu.uy>
 *
*/

#include <math.h>
/*#include <malloc.h>*/
#include <stdio.h>
#include "mex.h"
//#include "matrix.h"


/**
 * area_hist=areahist(label_image)
 */
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])

{
  size_t n,m;
  size_t mn;
  register size_t k;
  double *label_image;
  double *area_hist;
  /* temporary data */
  mxArray  *mx_label_areas;
  double *label_areas;

  size_t num_labels;
  
  /* 
   * check for correct # of input variables 
   */
  if (nrhs!=1){
    mexErrMsgTxt("There are exactly 1 input parameters allowed!");
    return;
  }

  /* 
   * get parameters
   */
  label_image = mxGetPr(prhs[0]);
  n = mxGetN(prhs[0]); 
  m = mxGetM(prhs[0]); 
  mn = n*m;

  /* count the labels */
  num_labels=0;
  for (k=0; k < mn; k++) {
    if (label_image[k]>num_labels)
      num_labels=label_image[k];
  }
  /* create label areas container (auto initialized to 0) */
  mx_label_areas = mxCreateDoubleMatrix(1,num_labels,mxREAL);
  label_areas = mxGetPr(mx_label_areas);

  /* compute label areas */
  for (k=0; k < mn; k++) {
    if (label_image[k] > 0.0) {
      label_areas[(int)label_image[k]-1] += 1.0;
    }
  }
  /* 
   * create output histogram and initialize it
   */
  plhs[0] = mxCreateDoubleMatrix(1,mn,mxREAL);
  area_hist = mxGetPr(plhs[0]);

  /* compute area histogram (discard 0 which is not a label) */
  for (k=0; k < num_labels; k++) {
    /* -1 because in Matlab index 0 is index 1 and must correspond to area 1*/
    area_hist[(int)label_areas[k]-1] += 1.0; 
  }

  /*
   * destroy temporary data
   */
  mxDestroyArray(mx_label_areas);
}
