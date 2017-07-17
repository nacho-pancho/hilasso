#include <mex.h>
#include "hilasso_admom.h"

#ifndef MATLAB_MEX_FILE
  #define MATLAB_MEX_FILE
#endif

/*
 * Input:
 * 0 x ............ (MxN) input vector or matrix of column vectors
 * 1 g ............ (Mx1) group indicators
 * 2 tau .......... (1x1) L1 penalty (lambda1) = tau
 * 3 fac12 ........ (1x1) L2 penalty (lambda2) = tau*lam1qlam2
 * 4 max_iter ..... (1x1) algorithm max. iterations (1000)
 * 5 tol .......... (1x1) algorithm tolerance (1e-8)
 * 6 c ............ (1x1) ADMOM parameter (50)
 * 
 * output:
 * 0 y ............ output
 */

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {

  /* pointers to raw data */
  /* input */
  double *x;
  double *g; 
  double lambda1(0),lambda2(0);
  double tol(1e-6), c(10);
  int max_iter(200);

  int* gi;
  int ng;

  /* output */
  double *y;

   /* auxiliary */
  register int i,j;
  int m,n;

   if (nrhs < 4) {
     mexErrMsgTxt("Short help goes here.\n");
       return;
   }

   /* input */
   x = mxGetPr(prhs[0]);   
   m = (int)mxGetM(prhs[0]);
   n = (int)mxGetN(prhs[0]);

   g = mxGetPr(prhs[1]);   
   int gm = (int)mxGetM(prhs[1]);
   int gn = (int)mxGetN(prhs[1]);
   if (gm == 1) {
     if (gn != m)
       mexErrMsgTxt("Bad dimension: g.\n");
   } else if (gn == 1) {
     if (gn != m)
       mexErrMsgTxt("Bad dimension: g.\n");
   } else {
       mexErrMsgTxt("Bad dimension: g.\n");
   }
   lambda1 = *mxGetPr(prhs[2]);    // tau
   lambda2 = (*mxGetPr(prhs[3]))*lambda1; // tau*
   if (nrhs > 4)
     max_iter = (int) *mxGetPr(prhs[4]);

   if (nrhs > 5)
     tol = *mxGetPr(prhs[5]);

   if (nrhs > 6)
     c = *mxGetPr(prhs[6]);
   
   //printf("nrhs=%d\tlambda1=%g\tlambda2=%g\tit=%d\ttol=%g\tc=%g\n",
	  //	  nrhs,lambda1,lambda2,max_iter,tol,c);

   /* 
    * initialize output 
    */
   plhs[0]=mxCreateDoubleMatrix(m,n,mxREAL);
   y=mxGetPr(plhs[0]);

   /*
    * find ranges of indices for each group
    * we assume groups are contiguous and ordered 
    */
   ng = int(g[m-1]);
   gi = new int[ng+1];
   j = 0;
   gi[j] = 0;
   for (i=0; i < m; i++)  {
     if (int(g[i]) > (j+1)) { // group has changed
       gi[++j] = i;
       //printf("group %d starts at %d\n",j+1,gi[j]);
     }
   }
   gi[ng]=m;
   //printf("\n");
   
   for ( j = 0 ; j < n ; j++ ) {
     for ( i = 0 ; i < ng ; i++ ) { 
       int mi = gi[i+1]-gi[i];
       // process each sample
       // L1 screening
       if (!scalar_soft_thres(&x[gi[i]],mi,lambda1,&y[gi[i]])) {
	 //putchar('s');
	 continue;
       }
       // L2 screening
       else if (!vector_soft_thres(&x[gi[i]],mi,lambda2,&y[gi[i]])) {
	 //putchar('v');
	 continue;
       }
       else { // proceed with ADMOM
	 //putchar('a');
	 bool res =
	   hilasso_admom(&x[gi[i]],mi,lambda1,lambda2,max_iter,tol,c,&y[gi[i]]);
	 /*
	   if (res)
	   putchar('.');
	   else
	   putchar('x');
	 */
       }
     }
   } 
   //putchar('\n');

   delete[] gi;
}
