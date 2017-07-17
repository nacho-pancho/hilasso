#include "mdls_linalg_types.h"

template<class T>
void mdls_deconstruct(const sc_matrix<T>& I, int w, int ov, 
		      const bool computeDC, const bool computeVAR,
		      sc_matrix<T>& X, 
		      sc_vector<T>& DC, 
		      sc_vector<T>& VAR )  // output 
{
   /* auxiliary */

  int m2,n2;

  int N;
  T dc,var;
  int d;
  
  register int i,j,i2,j2;
   
   int D = w*w; /* dimension of patches */
   int step = w - ov;
   int n = I.n;
   int m = I.m;

   /* these are temporary to compute the number of patches */
   j = (n-w)/step + 1; /* case when image breaks exactly */
   i = (m-w)/step + 1;
   /* 
    * if image doesn't break exactly into patches of this size
    * we pad it virtually.
    */
   if (((n-w) % step) != 0) {
     n2 = ((n-w)/step + 1)*step + w;
     j++;
   } else {
     n2 = n;
   }
   if (((m-w) % step) != 0) {
     m2 = ((m-w)/step + 1)*step + w;
     i++;
   } else {
     m2 = m;
   }   
   N = i*j; 
   //
   // allocate memory for returned data if not initialized yet
   // WARNING: we do not check that X and DC have the correct size!
   //
   if (X.n==0) {
     X.n = N;
     X.m = D;
     X.data = new T[D*N];
   }
   if (computeDC) {
     if (DC.n==0) {
       DC.data = new T[N];
       DC.n = N;
     }
   }
   if (computeVAR) {
     if (VAR.n==0) {
       VAR.data = new T[N];
       VAR.n = N;
     }
   }

   d = 0; /* since X is col-major, X 'wraps' to the next column automatically 
	     when d surpasses D */
   for ( i = 0 ; i <= (m2-w) ; i += step ) {
     for ( j = 0 ; j <= (n2-w) ; j += step ) {
       /* add a patch */
       dc = T(0);
       var = T(0);
       for ( j2 = 0; j2 < w ; j2++ ) {
          for ( i2 = 0 ; i2 < w ; i2++ ) {
	    T x = T();
	    if ( ( (i+i2) < m ) && ( (j+j2) < n ) )
	      x = I[m*(j+j2)+i+i2];
	    X[d++] = x;
	    dc  += x;
	    var += x*x;
	 }
       } 
       dc /= float(D);
       if (computeDC) {
	 DC[d/D-1] = dc;
       } // compute DC and VAR
       if (computeVAR) {
	 var = var/double(D) - dc*dc; // BIASED estimator
	 VAR[d/D-1] = var;	   
       } // computeVAR
     } // j
   }  // i
} // function
