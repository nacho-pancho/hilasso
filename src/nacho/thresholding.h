#include <cmath>
#include <cstring>

//
// soft thresholding. Return value is true if resulting L2 norm > 0
//
template<class T> bool scalar_soft_thres(const T* x, 
				       int n, 
				       const T& u, 
				       T* y) 
{
  register int i;
  bool a = false;
  for (i=0; i < n; i++) {
    if (x[i] > u) {
      y[i] = x[i] - u;
      a = true;
    } else if (x[i] < -u) {
      y[i] = x[i] + u;
      a = true;
    } else {
      y[i] = T(0);
    }
  } // main loop
  return a;
}

//
// vector thresholding. 
// Return value is true if resulting L2 norm > 0
//
template<class T> bool vector_soft_thres(const T* x, 
				       int n, 
				       const T& u, 
				       T* y) 
{
  // compute L2 norm
  register int i;
  T a(0);
  for (i=0; i < n; i++) {
    a += x[i]*x[i];
  }
  a = sqrt(a);
  if (a <= u) {
    // solution is zero
    memset(y,0,sizeof(T)*n);
    return false;
  } else {
    a = (a-u)/a;
    // solution is scaled down version of x
    for (i=0; i < n; i++) { 
      y[i] = a*x[i];
    }
    return true;
  }
}
