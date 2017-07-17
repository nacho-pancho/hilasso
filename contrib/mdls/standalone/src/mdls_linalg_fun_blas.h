#ifndef CSC_LINALG_FUN
#define CSC_LINALG_FUN
#include <ostream>

#include <gsl/gsl_cblas.h>

#ifdef DEBUG
#define my_assert(a) assert((a))
#else
#define my_assert(a) 
#endif

//
// basic operations
//
template<typename T> inline void sc_fill(const T& a, sc_vector<T>& u);
template<typename T> inline void sc_copy(const sc_vector<T>& u, sc_vector<T>& v);
template<typename T> inline void sc_scale(const T& a, sc_vector<T>& u);
template<typename T> inline void sc_scale(const T& a, const sc_vector<T>& u, sc_vector<T>& v);
template<typename T> inline void sc_add(const T& a, sc_vector<T>& u);
template<typename T> inline T    sc_sum(const sc_vector<T>& u);
template<typename T> inline int sc_amax(const sc_vector<T>& u);

//PEND template<typename T> inline void sc_add(const T& a, const sc_vector<T>& u, sc_vector<T>& v);

template<typename T> inline void sc_fill(const T& a, const sc_sparse_vector<T>& u);
template<typename T> inline void sc_copy(const sc_sparse_vector<T>& u, sc_sparse_vector<T>& v);
template<typename T> inline void sc_scale(const T& a, sc_sparse_vector<T>& u);
template<typename T> inline void sc_scale(const T& a, const sc_sparse_vector<T>& u, sc_sparse_vector<T>& v);
//PEND template<typename T> inline void sc_add(const T& a, sc_sparse_vector<T>& u);
//PEND template<typename T> inline void sc_add(const T& a, const sc_sparse_vector<T>& u, sc_sparse_vector<T>& v);

template<typename T> inline void sc_fill(const T& a, const sc_matrix<T>& A);
template<typename T> inline void sc_copy(const sc_matrix<T>& A, sc_matrix<T>& B);
template<typename T> inline void sc_scale(const T& a, sc_matrix<T>& A);
template<typename T> inline void sc_scale(const T& a, const sc_matrix<T>& A, sc_matrix<T>& B);
//PEND template<typename T> inline void sc_add(const T& a, sc_matrix<T>& A);
//PEND template<typename T> inline void sc_add(const T& a, const sc_matrix<T>& A, sc_matrix<T>& B);

template<typename T> inline void sc_fill(const T& a, const sc_sparse_matrix<T>& A); 
template<typename T> inline void sc_copy(const sc_sparse_matrix<T>& A, sc_sparse_matrix<T>& B); 
template<typename T> inline void sc_scale(const T& a, sc_sparse_matrix<T>& A);
template<typename T> inline void sc_scale(const T& a, const sc_sparse_matrix<T>& A, sc_sparse_matrix<T>& B);
//PEND template<typename T> inline void sc_add(const T& a, sc_sparse_matrix<T>& A);
//PEND template<typename T> inline void sc_add(const T& a, const sc_sparse_matrix<T>& A, sc_sparse_matrix<T>& B);
//
// conversion
//
template<typename T> void sc_full_to_sparse(const sc_vector<T>& vf, sc_sparse_vector<T>& vs);
template<typename T> void sc_sparse_to_full(const sc_sparse_vector<T>& vs, sc_vector<T>& vf);
//PEND template<typename T> void sc_pack(const sc_sparse_vector<T>& As);
//PEND template<typename T> void sc_full_to_sparse(const sc_matrix<T>& vf, sc_sparse_matrix<T>& vs);
//PEND template<typename T> void sc_sparse_to_full(const sc_sparse_matrix<T>& vs, sc_matrix<T>& vf);
//PEND template<typename T> void sc_pack(const sc_sparse_matrix<T>& As);

//
// addition
//
template<typename T>  void sc_add(const sc_vector<T>& u, const sc_vector<T>& v, sc_vector<T>& w);
template<typename T>  void sc_scaled_add(const T& a, const sc_vector<T>& u, const T& b, const sc_vector<T>& v, const sc_vector<T>& w);
//PEND template<typename T> inline T sc_add(const sc_vector<T>& u, const sc_sparse_vector<T>& v, const sc_vector<T>& w);
//PEND template<typename T> inline T sc_add(const sc_sparse_vector<T>& u, const sc_vector<T>& v, const sc_vector<T>& w);
template<typename T> void sc_add(const sc_sparse_vector<T>& u, const sc_sparse_vector<T>& v, sc_vector<T>& w);
template<typename T> void sc_scaled_add(const T& a, const sc_sparse_vector<T>& u, const T& b, const sc_sparse_vector<T>& v, sc_vector<T>& w);
template<typename T> void  sc_scaled_add(const T& a, const sc_vector<T>& v, sc_vector<T>& w);


//
// dot product
//
template<typename T> inline T sc_dot(const sc_vector<T>& u, const sc_vector<T>& v);
template<typename T> inline T sc_dot(const sc_vector<T>& u, const sc_sparse_vector<T>& v);
template<typename T> inline T sc_dot(const sc_sparse_vector<T>& u, const sc_vector<T>& v);
template<typename T> inline T sc_dot(const sc_sparse_vector<T>& u, const sc_sparse_vector<T>& v);

template<typename T> inline T sc_dot(const sc_matrix<T>& A, int irow, const sc_vector<T>& v);
template<typename T> inline T sc_dot(const sc_vector<T>& v, const sc_matrix<T>& A, int icol);
template<typename T> inline T sc_dot(const sc_matrix<T>& A, int irow, const sc_matrix<T>& B, int jcol);
//PEND template<typename T> inline T sc_dot(const sc_matrix<T>& A, int irow, const sc_sparse_vector<T>& v);
template<typename T> inline T sc_dot(const sc_sparse_vector<T>& v, const sc_matrix<T>& A, int icol);
//PEND template<typename T> inline T sc_dot(const sc_matrix<T>& A, int irow, const sc_matrix<T>& B, int jcol);

//
// matrix-vector product
//
template<typename T> inline void sc_mul(const sc_matrix<T>& A, const sc_vector<T>& u, sc_vector<T>& v);
template<typename T> inline void sc_mul(const sc_matrix<T>& A, bool At, const sc_vector<T>& u, sc_vector<T>& v);
template<typename T> inline void sc_mul(const sc_matrix<T>& A, const sc_sparse_vector<T>& u, sc_vector<T>& v);
//PENDING template<typename T> inline void sc_mul(const sc_sparse_matrix<T>& A, const sc_vector<T>& u, sc_vector<T>& v);
//PENDING template<typename T> inline void sc_mul(const sc_sparse_matrix<T>& A, const sc_sparse_vector<T>& u, sc_vector<T>& v);
//PENDING template<typename T> inline void sc_mul(const sc_sparse_matrix<T>& A, const sc_sparse_vector<T>& u, sc_sparse_vector<T>& v);
template<typename T> inline void sc_mul(const sc_vector<T>& u, const sc_matrix<T>& A, sc_vector<T>& v);
template<typename T> inline void sc_mul(const sc_sparse_vector<T>& u, const sc_matrix<T>& A, sc_vector<T>& v);
//PENDING template<typename T> inline void sc_mul(const sc_vector<T>& u, const sc_sparse_matrix<T>& A, sc_vector<T>& v);
//PENDING template<typename T> inline void sc_mul(const sc_sparse_vector<T>& u, const sc_sparse_matrix<T>& A, sc_vector<T>& v);
//PENDING template<typename T> inline void sc_mul(const sc_sparse_vector<T>& u, const sc_sparse_matrix<T>& A, sc_sparse_vector<T>& v);
//
// partial matrix-vector product: only compute rows which are already nonzero in target vector
//
template<typename T> inline void sc_mul_partial(const sc_matrix<T>& A, const sc_vector<T>& u, sc_sparse_vector<T>& v);
template<typename T> inline void sc_mul_partial(const sc_matrix<T>& A, const sc_sparse_vector<T>& u, sc_sparse_vector<T>& v);

//
// matrix-matrix
//
template<typename T> void sc_mul(const sc_matrix<T>& A, bool At, const sc_matrix<T>& B, bool Bt, sc_matrix<T>& C);
template<typename T> inline void sc_gram(const sc_matrix<T>& A, sc_matrix<T>& G);
template<typename T> inline void sc_cov(const sc_matrix<T>& A, sc_matrix<T>& S);
//
// display
//
template <typename T> void sc_dump(const sc_matrix<T>& A, std::ostream& out);
template <typename T> void sc_dump(const sc_vector<T>& v, std::ostream& out);
template <typename T> void sc_print(const sc_matrix<T>& A, const char* name);
template <typename T> void sc_print(const sc_vector<T>& v, const char* name);
//template <typename T> void sc_dump(const sc_matrix<T>& A, std::ostream& out);
//template <typename T> void sc_dump(const sc_vector<T>& v, std::ostream& out);
template <typename T> void sc_print(const sc_sparse_vector<T>& v, const char* name);
template <typename T> void sc_print(const sc_sparse_matrix<T>& A, const char* name);

//
// linear algebra basics
//
template<typename T> T sc_norm(const sc_vector<T>& u);
template<typename T>  T sc_norm(const sc_sparse_vector<T>& u);

template<typename T>  T sc_dist(const sc_vector<T>& u, const sc_vector<T>& v);
template<typename T>  T sc_dist(const sc_sparse_vector<T>& u, const sc_sparse_vector<T>& v);

//
// --------------------------------------------------------------------------------
//

//
// O(N) operations
//
//
// unary ops
//


template<typename T>
inline void sc_fill(const T& a, sc_vector<T>& u) { 
  register int i = 0;
  register int n = u.n;
  for ( ; i < n; i++) 
    u[i] = a;
}

template<typename T>
inline T sc_sum(const sc_vector<T>& u) { 
  register int i = 0;
  register int n = u.n;
  T a(0.0);
  for ( ; i < n; i++) 
    a += u[i];
  return a;
}

template<> inline int sc_amax<double>(const sc_vector<double>& u) {
  return cblas_idamax(u.n,u.data,1);
}

// copy
template<typename T>
inline void sc_copy(const sc_vector<T>& u, sc_vector<T>& v) { 
  register int i = 0;
  register int n = u.n;
  my_assert(u.n == v.n);
  for ( ; i < n; i++) 
    v[i] = u[i];
}

// destructive scaling
template<typename T>
inline void sc_scale(const T& a, sc_vector<T>& u) { 
  register int i = 0;
  register int n = u.n;
  for ( ; i < n; i++) 
    u[i] *= a;
}

// non-destructive scale
template<typename T>
inline void sc_scale(const T& a, const sc_vector<T>& u, sc_vector<T>& v) { 
  register int i = 0;
  register int n = u.n;
  my_assert(u.n == v.n);
  for ( ; i < n; i++) 
    v[i] = u[i]*a;
}

//
// sparse vector
//
template<typename T>
inline void sc_fill(const T& a, const sc_sparse_vector<T>& u) { 
  register int i = 0;
  register int n = u.nzmax;
  for ( ; i < n; i++) 
    u.values[i] = a;
}

// copy
template<typename T>
inline void sc_copy(const sc_sparse_vector<T>& u, sc_sparse_vector<T>& v) { 
  register int i = 0;
  register int n = u.nzmax;
  my_assert(u.n == v.n);
  my_assert(u.nzmax <= v.nzmax);
  for ( ; i < n; i++)  {
    v.indexes[i] = u.indexes[i];
    v.values[i] = u.values[i];
  }
}

template<typename T>
inline void sc_scale(const T& a, sc_sparse_vector<T>& u) { 
  register int i = 0;
  register int n = u.nzmax;
  for ( ; i < n; i++) 
    u.values[i] *= a;
}

template<typename T>
inline void sc_scale(const T& a, const sc_sparse_vector<T>& u, sc_sparse_vector<T>& v) { 
  register int i = 0;
  register int n = u.nzmax;
  my_assert(u.n == v.n);
  my_assert(u.nzmax <= v.nzmax);
  for ( ; i < n; i++)  {
    v.indexes[i] = u.indexes[i];
    v.values[i] = a*u.values[i];
  }
}


//
// matrix
//
template<typename T>
inline void sc_fill(const T& a, const sc_matrix<T>& A) { 
  register int i = 0;
  register int n = A.n * A.m;
  for ( ; i < n; i++) 
    A[i] = a;
}

// copy
template<typename T>
inline void sc_copy(const sc_matrix<T>& A, sc_matrix<T>& B) { 
  register int i = 0;
  my_assert(A.n == B.n);
  my_assert(A.m == B.m);
  register int n = A.n * A.m;
  for ( ; i < n; i++) 
    B[i] = A[i];
}

// destructive scaling
template<typename T>
inline void sc_scale(const T& a, sc_matrix<T>& A) { 
  register int i = 0;
  register int n = A.n *A.m;
  for ( ; i < n; i++) 
    A[i] *= a;
}

// non-destructive scale
template<typename T>
inline void sc_scale(const T& a, const sc_matrix<T>& A, sc_matrix<T>& B) { 
  register int i = 0;
  my_assert(A.n == B.n);
  my_assert(A.m == B.m);
  register int n = A.n * A.m;
  for ( ; i < n; i++) 
    B[i] = a*A[i];
}

//
// sparse matrix
//
template<typename T>
inline void sc_fill(const T& a, const sc_sparse_matrix<T>& A) { 
  register int i = 0;
  register int n = A.nzmax;
  for ( ; i < n; i++) 
    A.values[i] = a;
}

// copy
template<typename T>
inline void sc_copy(const sc_sparse_matrix<T>& A, sc_sparse_matrix<T>& B) { 
  register int i = 0;
  my_assert(A.n == B.n);
  my_assert(A.m == B.m);
  my_assert(A.nzmax <= B.nzmax);
  register int n = A.nzmax;
  for ( ; i < n; i++) {
    B.values[i] = A.values[i];
    B.row_indexes[i] = A.row_indexes[i];
  }
  n = A.n;
  for (i=0 ; i < n; i++) {
    B.col_offsets[i] = A.col_offsets[i];
  }
}

// destructive scaling
template<typename T>
inline void sc_scale(const T& a, sc_sparse_matrix<T>& A) { 
  register int i = 0;
  register int n = A.nzmax;
  for ( ; i < n; i++) 
    A.values[i] *= a;
}

// non-destructive scale
template<typename T>
inline T sc_scale(const T& a, const sc_sparse_matrix<T>& A, sc_sparse_matrix<T>& B) { 
  register int i = 0;
  my_assert(A.n == B.n);
  my_assert(A.m == B.m);
  my_assert(A.nzmax <= B.nzmax);
  register int n = A.nzmax;
  for ( ; i < n; i++) {
    B.values[i] = a*A.values[i];
    B.row_indexes[i] = A.row_indexes[i];
  }
  n = A.n;
  for (i=0 ; i < n; i++) {
    B.col_offsets[i] = A.col_offsets[i];
  }
}

//
// conversion sparse-full
//
// lightweight: space must be allocated already in vs: it only takes values from vf!
template<typename T> 
void sc_full_to_sparse(const sc_vector<T>& vf, sc_sparse_vector<T>& vs) 
{
  register int i = 0, is = 0;
  for (i = 0; i < vf.n; i++) {
    if (vf[i] != T(0)) {
      if (is >= vs.nzmax) {
	std::cerr << "Error converting to sparse: too many nonzero elements." << std::endl;
	break;
      }
      vs.values[is] = vf[i];
      vs.indexes[is++] = i+1;
    }
  }
  if (is < vs.nzmax)
    vs.indexes[is] = 0; 
}

template<typename T> 
void sc_sparse_to_full(const sc_sparse_vector<T>& vs, sc_vector<T>& vf) 
{
  sc_fill(T(0),vf);
  for (register int i = 0; i < vs.nzmax; i++) {
    if (vs.indexes[i] == 0)
      break;
    vf[vs.indexes[i]-1] = vs.values[i];
  }
}

//
// binary operations: dot product
//

//
// addition
//
template<typename T> 
void sc_add(const T& a, const sc_vector<T>& v)
{
  register int i = v.n -1;
  while (i >= 0) {
    v[i--] += a;
  }
}

template<typename T> 
void sc_scaled_add(const T& a, const sc_vector<T>& v, sc_vector<T>& w)
{
  my_assert(v.n == w.n);
  register int i = v.n -1;
  while (i >= 0) {
    w[i] = v[i] + a; i--;
  }
}

template<typename T> 
void sc_add(const sc_vector<T>& u, const sc_vector<T>& v, sc_vector<T>& w)
{
  my_assert(u.n == v.n);
  my_assert(v.n == w.n);
  register int i = u.n -1;
  while (i >= 0) {
    w[i] = v[i] + u[i]; i--;
  }
}

template<typename T> 
void sc_scaled_add(const T& a, const sc_vector<T>& u, const T& b, const sc_vector<T>& v, sc_vector<T>& w)
{
  my_assert(u.n == v.n);
  my_assert(v.n == w.n);
  register int i = u.n -1;
  while (i >= 0) {
    w[i] = a*u[i] + b*v[i]; i--;
  }
}

template<typename T> 
void sc_add(const sc_sparse_vector<T>& u, const sc_sparse_vector<T>& v, sc_sparse_vector<T>& w)
{
  my_assert(u.n == v.n);
  my_assert(v.n == w.n);
  T tmp[v.n];
  sc_vector<T> vtmp(tmp,v.n);
  register int i = u.n -1;
  register int j;
  sc_fill(T(),vtmp);
  while (i >= 0) {
    j = u.values[i]-1;
    if (j >= 0)
      tmp[j] += u.values[i];
    j = v.values[i]-1;    
    if (j >= 0)
      tmp[j] += v.values[i];
    i--;
  }
  sc_full_to_sparse(vtmp,w);
}

template<typename T>  T sc_norm(const sc_vector<T>& u)
{
  register int i = u.n;
  register T a = T();
  while (i >=0) {
    a += u[i]*u[i];
    i--;
  }
  return sqrt(a);
}

template<typename T>  T sc_norm(const sc_sparse_vector<T>& u)
{
  register int i = u.nzmax;
  register T a = T();
  while (i >=0) {
    if (u.indexes[i] > 0)
      a += u.values[i]*u.values[i];
    i--;
  }
  return sqrt(a);
}

template<typename T> T sc_dist(const sc_sparse_vector<T>& u, const sc_sparse_vector<T>& v)
{
  my_assert(u.n == v.n);
  register T a;
  int dnzmax = (v.nzmax*2 < v.n)? v.nzmax*2: v.n;
  int tmpi[dnzmax];
  T tmpv[dnzmax];
  sc_vector<int> tt(tmpi,dnzmax);
  sc_fill(0,tt);
  sc_vector<T> tt2(tmpv,dnzmax);
  sc_fill(T(),tt2);
  sc_sparse_vector<T> vtmp(tmpv,tmpi,v.n,dnzmax);
  sc_scaled_add(T(-1.0),u,T(1.0),v,vtmp);
  a = sc_norm(vtmp);
  return a;
}

template<typename T> T sc_dist(const sc_vector<T>& u, const sc_vector<T>& v)
{
  my_assert(u.n == v.n);
  T tmp[v.n];
  sc_vector<T> vtmp(tmp,v.n);
  sc_scaled_add(T(-1.0),u,T(1.0),v,vtmp);
  return sc_norm(vtmp);
}

template<typename T> void sc_scaled_add(const T& a, const sc_sparse_vector<T>& u, const T& b, const sc_sparse_vector<T>& v, sc_sparse_vector<T>& w)
{
  my_assert(u.n == v.n);
  my_assert(v.n == w.n);
  T tmp[v.n];
  sc_vector<T> vtmp(tmp,v.n);
  register int i = u.n -1;
  register int j;
  sc_fill(T(),vtmp);
  while (i >= 0) {
    j = u.values[i]-1;
    if (j >= 0)
      tmp[j] += a*u.values[i];
    j = v.values[i]-1;    
    if (j >= 0)
      tmp[j] += b*v.values[i];
    i--;
  }
  sc_full_to_sparse(vtmp,w);
}

//PEND template<typename T> inline T sc_add(const sc_vector<T>& u, const sc_sparse_vector<T>& v, const sc_vector<T>& w);
//PEND template<typename T> inline T sc_add(const sc_sparse_vector<T>& u, const sc_vector<T>& v, const sc_vector<T>& w);
//PEND template<typename T> inline T sc_add(const sc_sparse_vector<T>& u, const sc_sparse_vector<T>& v, const sc_vector<T>& w);

// dot product of two vectors
template<typename T>
inline T sc_dot(const sc_vector<T>& u, const sc_vector<T>& v) { 
  register int i = 0;
  register int n = u.n;
  register T a = T(0);
  my_assert(u.n == v.n);
  for ( ; i < n; i++) 
    a += u[i]*v[i];
  return a;
}

template<typename T> inline T sc_dot(const sc_vector<T>& u, const sc_sparse_vector<T>& v)
{
  register int i = 0, n = v.nzmax;
  register T a(0);
  my_assert(u.n == v.n);
  for ( ; i < n; i++) {
    if (v.indexes[i] < 1) 
      break;
    else
      a += v.values[i] * u[v.indexes[i]-1];
  }
  return a;
}

template<typename T> inline T sc_dot(const sc_sparse_vector<T>& u, const sc_vector<T>& v)
{
  return sc_dot(v,u);
}

template<typename T> inline T sc_dot(const sc_sparse_vector<T>& u, const sc_sparse_vector<T>& v)
{
  my_assert(u.n == v.n);
  register int i, j, iu = 0, iv = 0;
  int head_idx = 0; // which of the two index sets has the lead: 0->u, 1->v  
  int umax = u.nzmax;
  int vmax = v.nzmax;
  register T a(0);
  while ( (iu < umax) && (iv < vmax)  && (u.indexes[iu] > 0) && (v.indexes[iv] > 0) ) {
    if (head_idx == 0) { // u leads
      i = u.indexes[iu];
      j = v.indexes[iv];
      while (j < i)
	j = v.indexes[++iv];
      if (j == i)
	a += u.values[iu++] * v.values[iv++];
      else if (j > i) {
	head_idx = 1; // now leading is v
	iu++;
	continue;
      } else { 
	break;
      }
    } else { // v leads
      i = v.indexes[iv];
      j = u.indexes[iu];
      while (j < i)
	j = u.indexes[++iu];
      if (j == i)
	a += u.values[iv++] * v.values[iu++];
      else if (j > i) {
	head_idx = 0; // now leading is u
	iv++;
	continue;
      } else { 
	break;
      }
    } // end v leads
  } // end main loop
  return a;
}

// dot product of i-th row of a matrix and a vector
template<typename T>
inline T sc_dot(const sc_matrix<T>& A, int irow, const sc_vector<T>& v) { 
  register int n = A.n;
  register int m = A.m;
  register int iA = irow;
  register int iv = 0;
  register T a = T(0);
  my_assert(A.n == v.n);
  for ( ; iv < n; iv++, iA += m) // adding m to a col-major matrix jumps to next element in row
    a += A[iA]*v[iv];
  return a;
}

// dot product of a vector with the i-th column of a matrix
template<typename T>
inline T sc_dot(const sc_vector<T>& v, const sc_matrix<T>& A, int icol) { 
  return sc_dot(v,A.column(icol));
}

template<typename T> inline T sc_dot(const sc_sparse_vector<T>& v, const sc_matrix<T>& A, int icol)
{
  return sc_dot(v,A.column(icol));
}

// dot product of the i-th row of a matrix  with the j-th column of another matrix
template<typename T>
inline T sc_dot(const sc_matrix<T>& A, int irow, const sc_matrix<T>& B, int jcol) { 
  register int Am = A.m;
  register int iA = irow;
  register int iB = B.m*jcol;
  register int fiB = iB+B.m;
  register T a = T(0);
  my_assert(A.n == B.m);
  for ( ; iB < fiB ; iB++, iA += Am) 
    a += A[iA]*B[iB];
  return a;
}

//
// matrix-vector products
//

// matrix-vector product 
template<typename T>
inline void sc_mul(const sc_matrix<double>& A, const sc_vector<double>& u, sc_vector<T>& v) { 
  sc_mul(A,false,u,v);
}

template<>
inline void sc_mul<double>(const sc_matrix<double>& A, bool At, const sc_vector<double>& u, sc_vector<double>& v) {    
  cblas_dgemv(CblasColMajor, // storage order
	      At? CblasTrans: CblasNoTrans, // transpose
	      At? A.n: A.m, // M
	      At? A.m: A.n, // N
	      1.0, // alpha 
	      A.data, // 
	      A.m, // lda
	      u.data,
	      1, // incX
	      0.0, // beta
	      v.data,
	      1); // incY
	      
}

template<>
inline void sc_mul<float>(const sc_matrix<float>& A, bool At, const sc_vector<float>& u, sc_vector<float>& v) {    
 cblas_sgemv(CblasColMajor, // storage order
	      At? CblasTrans: CblasNoTrans, // transpose
	      At? A.n: A.m, // M
	      At? A.m: A.n, // N
	      1.0f, // alpha 
	      A.data, // 
	      A.m, // lda
	      u.data,
	      1, // incX
	      0.0f, // beta
	      v.data,
	      1); // incY
	      
}

template<typename T>
inline void sc_mul(const sc_matrix<T>& A, const sc_sparse_vector<T>& u, sc_vector<T>& v) { 
  my_assert(A.n == u.n);
  my_assert(A.m == v.n);
  register int m = A.m;
  //register int n = A.n;
  register int i,j,iA;
  register T a;
  register int nzi = 0;
  sc_fill(T(0),v);
  for (nzi=0; nzi < u.nzmax; nzi++) {
    if (u.indexes[nzi] < 1)
      break;
    j = u.indexes[nzi]-1;
    a = u.values[nzi];
    iA = m*j;
    for (i=0; i < m; i++) {
      v[i] += a*A[iA+i];
    }
  }
}

template<typename T>
inline void sc_mul(const sc_vector<T>& u, const sc_matrix<T>& A, sc_vector<T>& v) { 
  sc_mul(A,true,u,v);
}

template<typename T>
inline void sc_mul(const sc_sparse_vector<T>& u, const sc_matrix<T>& A, sc_vector<T>& v) { 
  my_assert(A.m == u.n);
  my_assert(A.n == v.n);
  register int i,j,iA;
  register T a;
  register int m = A.m;
  register int n = A.n;
  register int nzi = 0;
  iA = 0;
  for (j = 0; j < n; j++) {
    a = T(0);
    for (nzi=0; nzi < u.nzmax; nzi++) {
      i =  u.indexes[nzi] - 1;
      a += A[iA+i] * u.values[nzi];
    }
    v[j]=a;
    iA += m;
  }
}

// partial matrix-vector product 
template<typename T>
inline void sc_mul_partial(const sc_matrix<T>& A, const sc_vector<T>& u, sc_sparse_vector<T>& v) { 
  my_assert(A.n == u.n);
  my_assert(A.m == v.n);
  register int n = A.n;
  register int i,j,k,iA;
  register T a;
  iA = 0;
  sc_fill(T(0),v); // support is kept in v.indexes, not touched by fill
  for (j=0; j < n; j++) {
    a = u[j];
    iA = A.m*j;
    for (i=0; i < v.nzmax; i++) {
      k = v.indexes[i]-1;
      v.values[i] += a*A[iA+k];
    }
  }
}
/*
template<typename T>
inline void sc_mul(const sc_matrix<T>& A, const sc_sparse_vector<T>& u, sc_sparse_vector<T>& v) { 
  my_assert(A.n == u.n);
  my_assert(A.m == v.n);
  register int m = A.m;
  //register int n = A.n;
  register int i,j,iA;
  register T a;
  register int nzi = 0;
  sc_fill(T(0),v);
  for (nzi=0; nzi < u.nzmax; nzi++) {
    j = u.indexes[nzi]-1;
    a = u.values[nzi];
    iA = m*j;
    for (i=0; i < m; i++) {
      v[i] += a*A[iA+i];
    }
  }
}
*/

//
// matrix-matrix product 
//
template<typename T>
void sc_mul(const sc_matrix<T>& A, bool At, const sc_matrix<T>& B, bool Bt, sc_matrix<T>& C) { 
  //
  //
  //
  int Am,An,Bm,Bn, dAi,dAj,dBi,dBj;
  int iA,iB,iC;
  register int i,j,k,m,n;
  register double a;
  int iA0, iB0;
  if (At) {
    Am = A.n;
    An = A.m;
    dAi = A.m;
    dAj = 1;
  } else {
    Am = A.m;
    An = A.n;
    dAi = 1;
    dAj = A.m;
  }
  if (Bt) {
    Bm = B.n;
    Bn = B.m;
    dBi = B.m;
    dBj = 1;
  } else {
    Bm = B.m;
    Bn = B.n;
    dBi = 1;
    dBj = B.m;
  }  
  m = C.m;
  n = C.n;
  my_assert(An == Bm);
  my_assert(Am == C.m);
  my_assert(Bn == C.n);
  iC = 0;
  iB0 = 0;
  for (j = 0; j < n; j++) {
    iA0 = 0;
    for (i = 0; i < m; i++) {
      a = T(0);
      iA = iA0;
      iB = iB0;
      for (k = 0; k < An; k++, iA += dAj, iB += dBi) {
	//std::cout << "iA=" << iA << " iB=" << iB << " iC=" << iC << std::endl;
	a += A[iA]*B[iB];
      }
      C[iC++]=a;
      iA0 += dAi;
    }
    iB0 += dBj;
  }
}

// gram matrix product: G=At*A
template<typename T>
inline void sc_gram(const sc_matrix<T>& A, sc_matrix<T>& G) { 
  sc_mul(A,true,A,false,G);
}

// convariance-like matrix product : S=A*At
template<typename T>
inline void sc_cov(const sc_matrix<T>& A, sc_matrix<T>& S) { 
  sc_mul(A,false,A,true,S);
}


template <typename T> void sc_dump(const sc_matrix<T>& A, std::ostream& out) {
   out.precision(10);
   int _n = A.n;
   int _m = A.m;
   for (int i = 0; i<_m; ++i) {
      for (int j = 0; j<_n; ++j) {
         out << A[j*_m+i] << " ";
      }
      out << std::endl;
   }
   out << std::endl;
}

template <typename T> void sc_print(const sc_matrix<T>& A, const char* name) {
   std::cerr << name << std::endl;
   std::cerr.precision(5);
   std::cerr << std::fixed;
   int _n = A.n;
   int _m = A.m;
   std::cerr << _m << " x " << _n << std::endl;
   for (int i = 0; i<_m; ++i) {
      for (int j = 0; j<_n; ++j) {
         std::cerr << A[j*_m+i] << " ";
      }
      std::cerr << std::endl;
   }
   std::cerr << std::endl;
}


template <typename T> void sc_dump(const sc_vector<T>& v, std::ostream& out)  {
   out.precision(10);
   int _n = v.n;
   for (int i = 0; i<_n; ++i) {
      out << v[i] << ' ';
   }
   out << std::endl;
};

template <typename T> void sc_print(const sc_vector<T>& v, const char* name)  {
   std::cerr.precision(5);
   std::cerr << std::fixed;
   std::cerr << name << std::endl;
   int _n = v.n;
   for (int i = 0; i<_n; ++i) {
      std::cerr << v[i] << " ";
   }
   std::cerr << std::endl;
};

template <typename T> void sc_print(const sc_sparse_vector<T>& v, const char* name)  {
   std::cerr.precision(5);
   std::cerr << std::fixed;
   std::cerr << name << " (sparse)" << std::endl;
   int i,nzi,j;
   for (i = 0,nzi = 0; i < v.n; ++i) {
     if (nzi >= v.nzmax) {
       std::cerr << T(0) << " ";
       continue;
     }
     j = v.indexes[nzi];
     if (j < 1) {
       std::cerr << T(0) << " ";
       nzi = v.nzmax + 1; // from now on all are zeros
       continue;
     }
     if (i == (j-1)) {
       std::cerr << v.values[nzi] << " ";
       nzi++;
     } else {
       std::cerr << T(0) << " ";
     }
   }
   std::cerr << std::endl;
};

template <typename T> void sc_print(const sc_sparse_matrix<T>& A, const char* name) {
   std::cerr << name << std::endl;
   std::cerr.precision(5);
   std::cerr << std::fixed;
   int _n = A.n;
   int _m = A.m;   
   int nzi;
   std::cerr << _m << " x " << _n << " (sparse)" << std::endl;
   for (int i = 0; i<_m; ++i) {
      for (int j = 0; j<_n; ++j) {
	for (nzi = A.col_offsets[j]; nzi < A.col_offsets[j+1]; nzi++) {
	  if ((A.row_indexes[nzi]-1) == i) {
	    std::cerr << A.values[nzi] << " ";
	    break;
	  }
	}
        if (nzi == A.col_offsets[j+1]) 
	    std::cerr << T(0) << " ";
      }
      std::cerr << std::endl;
   }
   std::cerr << std::endl;
};

#undef my_assert
#endif
