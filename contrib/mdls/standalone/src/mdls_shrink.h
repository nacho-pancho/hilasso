#ifndef MDLS_SHRINK
#define MDLS_SHRINK
#include "mdls_linalg.h"

// destructive
template <typename T> inline void sc_soft_threshold(const T threshold, sc_vector<T>& u);
template <typename T> inline void sc_soft_threshold(const T threshold, sc_sparse_vector<T>& u);
// non-destructive
template <typename T> inline void sc_soft_threshold(const T threshold, const sc_vector<T>& u, sc_vector<T>& v);
template <typename T> inline void sc_soft_threshold(const T threshold, const sc_vector<T>& u, sc_sparse_vector<T>& v);
// partial update only on the already non-zero elements of the target vector
template <typename T> inline void sc_soft_threshold_partial(const T threshold, const sc_vector<T>& u, sc_sparse_vector<T>& v);

template <typename T> 
inline void sc_soft_threshold(const T threshold, sc_vector<T>& u) {
  register int n = u.n;
  for (register int i = 0; i < n ; ++i) {
    if (u[i] > threshold) {
      u[i] -= threshold;
    } else if (u[i] < -threshold) {
      u[i] += threshold;
    } else {
      u[i] = T();
    }
  }
}

template <typename T> 
inline void sc_soft_threshold(const sc_vector<T>& tv, sc_vector<T>& u) {
  register int n = u.n;
  assert(tv.n == u.n);
  for (register int i = 0; i < n ; ++i) {
    if (u[i] > tv[i]) {
      u[i] -= tv[i];
    } else if (u[i] < -tv[i]) {
      u[i] += tv[i];
    } else {
      u[i] = T();
    }
  }
}

template <typename T> 
inline void sc_soft_threshold(const T threshold, sc_sparse_vector<T>& u) {
  register int n = u.nmax;
  for (register int i = 0; i < n ; ++i) {
    if (u.values[i] > threshold) {
      u.values[i] -= threshold;
    } else if (u.values[i] < -threshold) {
      u.values[i] += threshold;
    } else {
      u.values[i] = T();
    }
  }
}

template <typename T> 
inline void sc_soft_threshold(const T threshold, const sc_vector<T>& u, sc_vector<T>& v) {
  register int n = v.n;
  assert(u.n == v.n);
  for (register int i = 0; i < n ; ++i) {
    if (u[i] > threshold) {
      v[i] = u[i] - threshold;
    } else if (u[i] < -threshold) {
      v[i] = u[i] + threshold;
    } else {
      v[i] = T();
    }
  }
}

template <typename T> 
inline void sc_soft_threshold(const T threshold, const sc_vector<T>& u, sc_sparse_vector<T>& v) {
  assert(u.n == v.n);
  register int n = u.n;
  register int i =0, j=0;

  for (; i < n ; ++i) {
    if (u[i] > threshold) {
      assert(j < v.nzmax);
      v.values[j] = u[i] - threshold;
      v.indexes[j++]=i+1;
    } else if (u[i] < -threshold) {
      assert(j < v.nzmax);
      v.values[j] = u[i] + threshold;
      v.indexes[j++]=i+1;
    }
  }
}

template <typename T> 
inline void sc_partial_soft_threshold(const T threshold, const sc_vector<T>& u, sc_sparse_vector<T>& v) {
  assert(u.n == v.n);
  register int n = v.nzmax;
  register int i =0,j=0;

  for (; i < n ; ++i) {
    if (v.indexes[i] < 1)
      continue;
    j = v.indexes[i]-1;
    if (u[j] > threshold) {
      v.values[i] = u[j] - threshold;
    } else if (u[i] < -threshold) {
      v.values[i] = u[j] + threshold;
    } else {
      v.values[i] = T();
    }
  }
}

#endif
