#include <math.h>

#define MACRO_clip_coords(x,n) (  ( (((x) >=0) ? (x) : 0 ) < (n) )? (x) : ((n)-1)  )

void mdls_get_patches(const T* I, int m, int n, const int* grid, int N, // input
		      T* X, T* DC, T* VAR) // output
{
  double dc,var;
  int d;
  bool computeDC,computeVAR;
  computeDC = (DC != NULL);
  computeVAR = (VAR != NULL);

   int D  = w*w; /* dimension of patches */
   int w2 = w/2; /* patch radius rounded to nearest smallest integer */
   int d  = 0; /* since X is col-major, X 'wraps' to the next column automatically 
	     when d surpasses D */
   for (int k = 0; k < (2*N); ) {
     int j = grid[k++]-w2; // j indexes the 'x' direction
     int i = grid[k++]-w2; // i indexes the 'y' dir
     /* add a patch */
     dc = 0.0;
     var = 0.0;
     for (int j2 = 0; j2 < w ; j2++ ) {
       for (int i2 = 0 ; i2 < w ; i2++ ) {
	 T x;
	 int i3 = MACRO_clip_coords(i+i2,m);
	 int j3 = MACRO_clip_coords(j+j2,n);
	 X[d++] = x = I[ m*j3 + i3 ];	 
	 dc += x;
	 var += x*x;
       } // i2
     } // j2, inner loop to copy patch data
     if (computeDC) {
       /* compute DC */
       dc /= double(D);
       DC[d/D-1] = dc;
       if (computeVAR) {
	 var = var/double(D) - dc*dc; // BIASED estimator
	 VAR[d/D-1] = var;
       }
     }
   } // for each  grid location
}
