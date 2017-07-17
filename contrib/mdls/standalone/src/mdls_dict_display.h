#ifndef MDLS_DICT_DISPLAY
#define MDLS_DICT_DISPLAY
#include <cmath>

template<typename T, typename P>
  void mdls_dict_display(const sc_matrix<T>& D, sc_matrix<P>& I,
			 int margin=1, P bg=0) 
{
  int K = D.n;
  int w = int(sqrt(double(D.m)));
  int ng = (int) ceil(sqrt(double(K))); // atom grid width
  int mg = (int) ceil(double(K)/double(ng));    // atom grid height
  int ni = ng*w + margin*(ng-1);        // output image width
  int mi = mg*w + margin*(mg-1);        // output image height

/* #if DEBUG */
/*   std::cout << "K=" << K << " M=" << M << " w=" << w  */
/* 	    << " mg=" << mg << " ng=" << ng << " mi=" << mi << " ni=" << ni  */
/* 	    << std::endl; */
/* #endif */
    
  //
  // eventually, allocate space for output
  //
  if ((I.m != mi) || (I.n != ni)) {
    if (I.data != NULL) {
      I.free();
    }
    I.allocate(mi,ni);
  }
  sc_fill( bg, I);
  // 
  // fill
  //
  int k = 0;
  for (int ig = 0; ig < mg; ig++) {
    for (int jg = 0; (jg < ng) && (k < K); jg++, k++) {      
      int ii0=ig*(w+margin);
      int ii1=ii0+w;
      int ji0=jg*(w+margin);
      int ji1=ji0+w;
      // inner 2D loop
      sc_vector<T> Dk = D.column(k);
      int l = 0;
      for (int ji=ji0; ji <  ji1; ji++) {
	for (int ii=ii0; ii <  ii1; ii++) {
	  I(ii,ji) = (P) Dk[l++];
	}
      } // inner: paint atom
    }
  } // outer: for each atom
  return;
} // function

#endif
