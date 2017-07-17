#include <cmath>

/*
 * Input:
 * 0 X ............ patches (D x N)
 * 1 w ............ width of patches
 * 2 N ............ number of patches
 * 1 m ............ height of output image
 * 2 n ............ width of output image
 * 3 ov ........... overlap of patches
 % 4 DC ........... optional DC, (1 x N)
 *
 * output:
 * 0 I ............ Output image.
 * 0 R ............ Working space (same size as img)
 */
template<typename T>
void mdls_reconstruct(const T* X, int w, int N,
		     int m, int n, int ov, const T* DC, 
		     T*& I, T* R) {
   /* auxiliary */
   bool internalR = (R == NULL);

   /* dimensions */
   int D  = w*w;
   int L  = m*n;
   int step = w - ov;

   /* initialize auxiliary normalization matrix */
   if (internalR) {
     R = new T[L];
   }
   for (int i = 0; i < L ; i++) {
     R[i] = 0.0; I[i] = 0.0;
   }
   int d = 0;
   double dc;
   for (int i = 0 ; i <= (m-w) ; i += step ) {
     for (int j = 0 ; j <= (n-w) ; j += step ) {
       /* add a patch */
       dc = (DC!=NULL) ? DC[d/D] : 0.0;
       for (int j2 = 0; j2 < w ; j2++ ) {
          for (int i2 = 0 ; i2 < w ; i2++ ) {
	   int idx = m*(j+j2)+i+i2;
	   I[idx] += (X[d++] + dc); 
	   R[idx] += 1.0 ;
	 }
       } 
     }
   } 

   /* normalize */
   for (int i = 0; i < L ; i++) {
     if (R[i] > 0.0)
       I[i] /= R[i];
   }
   
  /* destroy temporary space */
   if (internalR) {
     delete[] R;
   }
}

template<typename T>
void mdls_reconstruct_with_weight(const T* X, int w, int N, const T* W, int m, int n, int ov, const T* DC, T* I, T* R) {
   /* auxiliary */
   bool internalR = (R == NULL);

   /* dimensions */
   int L = m*n;
   int D = w*w;
   int step = w - ov;

   /* initialize normalization matrix */
   if (internalR) {
     R = new T[L];
   }
   for (int i = 0; i < L ; i++) {
     R[i] = 0.0; I[i] = 0.0;
   }
   int d = 0;
   double dc;
   for (int i = 0 ; i <= (m-w) ; i += step ) {
     for (int j = 0 ; j <= (n-w) ; j += step ) {
       /* add a patch */
       dc = (DC!=NULL) ? DC[d/D] : 0.0;
       int wd = 0;
       for (int j2 = 0; j2 < w ; j2++ ) {
          for (int i2 = 0 ; i2 < w ; i2++ ) {
	   int idx = m*(j+j2)+i+i2;
	   I[idx] += (X[d++] + dc) * W[wd]; 
	   R[idx] += W[wd++];
	 }
       } 
     }
   } 

   /* normalize */
   for (int i = 0; i < L ; i++) {
     if (R[i] > 0.0)
       I[i] /= R[i];
   }
   
  /* destroy temporary space */
   if (internalR) {
     delete[] R;
   }
}
