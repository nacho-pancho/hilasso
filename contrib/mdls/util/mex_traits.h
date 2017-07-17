#ifndef MEX_TRAITS_H
#define MEX_TRAITS_H
#include <mex.h>

template<class T>
struct mex_traits {
  typedef T prec_t;
  inline static mxClassID mx_class() { return mxClassID(0); }
};

template<> struct mex_traits<char> {
  typedef char prec_t;
  typedef short op_t;
  inline static mxClassID mx_class() { return mxCHAR_CLASS; }
};

template<> struct mex_traits<unsigned char> {
  typedef char prec_t;
  typedef short op_t;
  inline static mxClassID mx_class() { return mxUINT8_CLASS; }
};

template<> struct mex_traits<short> {
  typedef char prec_t;
  typedef short op_t;
  inline static mxClassID mx_class() { return mxINT16_CLASS; }
};

template<> struct mex_traits<unsigned short> {
  typedef unsigned short prec_t;
  typedef short op_t;
  inline static mxClassID mx_class() { return mxUINT16_CLASS; }
};

template<> struct mex_traits<int> {
  typedef int prec_t;
  typedef int op_t;
  inline static mxClassID mx_class() { return mxINT32_CLASS; }
};

template<> struct mex_traits<unsigned int> {
  typedef unsigned int prec_t;
  typedef int op_t;
  inline static mxClassID mx_class() { return mxUINT32_CLASS; }
};

template<> struct mex_traits<long> {
  typedef long prec_t;
  typedef long op_t;
  inline static mxClassID mx_class() { return mxINT64_CLASS; }
};

template<> struct mex_traits<unsigned long> {
  typedef long prec_t;
  typedef long op_t;
  inline static mxClassID mx_class() { return mxUINT64_CLASS; }
};

template<> struct mex_traits<float> {
  typedef float prec_t;
  typedef double op_t;
  inline static mxClassID mx_class() { return mxSINGLE_CLASS; }
};

template<> struct mex_traits<double> {
  typedef double prec_t;
  typedef double op_t;
  inline static mxClassID mx_class() { return mxDOUBLE_CLASS; }
};

#endif
