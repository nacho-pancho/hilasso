#ifndef CSC_LINALG_TYPES
#define CSC_LINALG_TYPES
#include <cassert>
#include <iostream>
#include <iomanip>

#ifdef DEBUG
#define my_assert(a) assert((a))
#else
#define my_assert(a) 
#endif

template<typename T>
struct sc_vector {
  //
  // data
  //
  int n;
  T* data;  
  typedef T type;
  //
  // operations
  //
  sc_vector(): n(0), data(NULL) {}
  sc_vector(T* _data, int _n): n(_n),data(_data) {}
  inline void allocate(int _n) {
    n = _n; data = new T[n];
  }
  inline void free() { 
    if (data != NULL) {
      delete[] data;
      data = NULL;
    }
    n = 0; 
  }
  inline T& operator[](int i) { my_assert(i < n); return data[i]; }
  inline T operator[](int i) const { my_assert(i < n); return data[i]; }
};

/**
 * sparse vector implementation
 * vectors are of dimension n
 * data is stored in two arrays of size nzmax <= n
 * - indexes
 * - values
 * 
 * a dimension n, a max. number of nonzero elements nzmax
 * a values buffer of size nzmax and an indexes buffer of size nzmax.
 * 
 * nonzero 
 */
template<typename T>
struct sc_sparse_vector 
{
  //
  // data
  //
  int n, nzmax;
  T* values;
  int* indexes; // 1-offset
  //
  // operations
  //
  sc_sparse_vector(): 
  n(0), nzmax(0), values(NULL), indexes(NULL) {}

sc_sparse_vector(T* _values, int* _indexes, int _n, int _nzmax): 
  n(_n), nzmax(_nzmax), values(_values), indexes(_indexes) {}


  inline void allocate(int _n, int _nzmax) {
    int i;
    n = _n; nzmax = _nzmax;     
    values = new T[nzmax]; indexes = new int[nzmax];
    for (i=0; i < nzmax; i++) {
      indexes[i] = -1;
    }
  }
  //
  // index access: SLOW
  //
  T operator[](int i) const {
    for (int ii=0; (indexes[ii] >= 0) && (ii < nzmax); ii++) {
      if (indexes[ii]==i) {
	return values[ii];
      }
    }
    return T(0);
  }

  void free() {
    if (values != NULL) {
      delete[] values;
      values = NULL;
      delete[] indexes;
      indexes = NULL;
    }
    n = 0; nzmax = 0;
  }

  int nnz() const {
    int ii;
    for (ii=0; (indexes[ii] >= 0) && (ii < nzmax); ii++) {}
    return ii;
  }

  // these are more heavy now, better use value and index functions for 'smart' usage
  //T& operator[](int i); 
  //const T& operator[](int i) const; 

  inline T& val(int i) { my_assert(i < nzmax); values[i]; }
  inline const T& val(int i) const { my_assert(i < nzmax); values[i]; }
  inline T& idx(int i) { my_assert(i < nzmax); idx[i]; }
  inline const T& idx(int i) const { my_assert(i < nzmax); idx[i]; }

};

template<typename T>
struct sc_matrix {
  //
  // data
  //
  int m,n;
  T* data;

  //
  // operations
  //
 sc_matrix(): m(0),n(0),data(NULL) {}
 sc_matrix(T* _data,int _m,int _n): m(_m),n(_n),data(_data) {}
  inline T& operator[](int i) { my_assert(i < (m*n)); return data[i]; }
  inline const T& operator[](int i) const { my_assert(i < (m*n)); return data[i]; }  
  inline T& operator()(int i,int j) const { my_assert((i < m)&&(j < n)); return data[m*j+i]; }  
  inline sc_vector<T> column(int i) { my_assert(i < n); return sc_vector<T>(&data[m*i],m); }
  inline const sc_vector<T> column(int i) const { my_assert(i < n); return sc_vector<T>(&data[m*i],m); }
  inline void allocate(int _m, int _n) {
    m = _m; n = _n; data = new T[m*n];
  }
  void free() {
    if (data != NULL) {
      delete[] data;
      data = NULL;
    }
    n = 0; m = 0;
  }

};

/**
 * similar to Matlab format, but the number of nonzeros are imposed
 * on a per column basis, which makes things easier.
 * there are two matrices of size nzmax x n:
 *
 * row index matrix
 * values matrix
 *
 */
template<typename T>
struct sc_sparse_matrix {
  //
  // data
  //
  int m,n;
  T* values;
  int* row_indexes;
  int col_nzmax;
  int nzmax; // total number of nonzero elements
  //
  // operations
  //
  sc_sparse_matrix(): 
    m(0),n(0),values(NULL),row_indexes(NULL),col_nzmax(0) {}

  sc_sparse_matrix(T* _values, 
		   int* _row_indexes, 
		   int _m,
		   int _n, 
		   int _col_nzmax): 
  m(_m),n(_n),values(_values),row_indexes(_row_indexes), col_nzmax(_col_nzmax) {}

  inline void allocate(int _m, int _n, int _col_nzmax) {
    m = _m; 
    n = _n; 
    col_nzmax   = _col_nzmax; 
    nzmax       = n*col_nzmax;
    values      = new T[col_nzmax*n]; 
    row_indexes = new int[col_nzmax*n];
  }

  inline void free() {
    if (values != NULL) { 
      delete[] values; 
      values = NULL;
    }
    if (row_indexes != NULL) { 
      delete[] row_indexes; 
      row_indexes = NULL;
    }
    col_nzmax = nzmax = n = m = 0;
  }

  inline sc_sparse_vector<T> column(int i) { 
    my_assert(i < n); 
    return sc_sparse_vector<T>(&values[col_nzmax*i], 
			       &row_indexes[col_nzmax*i], 
			       m, 
			       col_nzmax); 
  }

  inline const sc_sparse_vector<T> column(int i) const { 
    my_assert(i < n); 
    return sc_sparse_vector<T>(&values[col_nzmax*i], 
			       &row_indexes[col_nzmax*i], 
			       m, 
			       col_nzmax); 
  }

  inline const T operator()(int i, int j) const {
    const sc_sparse_vector<T> v = this->column(j);
    return v[i];
  }

};
#undef my_assert
#endif
