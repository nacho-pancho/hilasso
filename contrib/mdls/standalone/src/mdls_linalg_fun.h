#ifndef CSC_LINALG_FUN
#define CSC_LINALG_FUN
#include <ostream>
#include <cstdlib>

#ifdef DEBUG
#define my_assert(a) assert((a))
#else
//#define my_assert(a) 
#define my_assert(a) assert((a))
#endif

//
// --------------------------------------------------------------------------------
//

template<typename T> inline void sc_rank_one_update(const T& alpha, const sc_vector<T>& u, const sc_vector<T>& v, sc_matrix<T>& A)
{
  my_assert(u.n == A.m);
  my_assert(v.n == A.n);
  register int k = 0;
  register int N = v.n;
  register int M = u.n;
  for (int j=0; j < N; j++) {
    T avj = alpha*v.data[j];
    for (int i=0; i < M; i++) {
      A.data[k++] += avj*u.data[i];
    }
  }
}

template<typename T> inline void sc_rank_one_update(const T& alpha, const sc_vector<T>& u, const sc_sparse_vector<T>& v, sc_matrix<T>& A)
{
  my_assert(u.n == A.m);
  my_assert(v.n == A.n);
  int nnz = v.nzmax;
  int j;
  for (int k=0; (k < nnz) && ((j=v.indexes[k])>=0); k++) {
    sc_vector<T> Aj = A.column(j);
    sc_scaled_add(alpha*v.values[k],u,Aj);
  }
}


template<typename T>
inline bool sc_isnan(const sc_vector<T>& u) { 
  register int i = 0;
  register int n = u.n;
  for ( ; i < n; i++) 
    if (isnan(u.data[i]) )
      return true;
  return false;
}

template<typename T>
inline bool sc_isnan(const sc_matrix<T>& A) { 
  register int i = 0;
  register int n = A.n * A.m;
  for ( ; i < n; i++) 
    if (isnan(A.data[i]) )
      return true;
  return false;
}

template<typename T>
inline bool sc_isnan(const sc_sparse_vector<T>& u) { 
  register int i;
  register int n = u.nzmax;
  for (i=0 ;  (u.indexes[i] >= 0) && (i < n); i++) 
    if ( isnan(u.values[i]) )
      return true;
  return false;
}

template<typename T>
inline bool sc_isnan(const sc_sparse_matrix<T>& A) { 
  register int j = 0;
  register int n = A.n;
  for ( ; j < n; j++) {
    sc_sparse_vector<T> Aj = A.column(j);
    if (isnan(Aj))
      return true;
  }
  return false;
}

template<typename T>
inline void sc_fill(const T& a, sc_vector<T>& u) { 
  register int i = 0;
  register int n = u.n;
  for ( ; i < n; i++) 
    u[i] = a;
}

template<typename T>
inline void sc_clear(const T& a, sc_vector<T>& u) { 
  sc_fill(T(0),u);
}

template<typename T>
inline void sc_clear(const T& a, sc_matrix<T>& A) { 
  sc_fill(T(0),A);
}

template<typename T>
inline void sc_clear(const T& a, sc_sparse_vector<T>& u) { 
  if(u.nzmax>0) {
    u.indexes[0]= -1;
  }
}

template<typename T>
inline void sc_clear(const T& a, sc_sparse_matrix<T>& A) { 
  for (int j=0; j < A.n; j++) {
    sc_sparse_vector<T> Aj = A.column(j);
    sc_clear(Aj);
  }
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

template<typename T> inline int sc_amax(const sc_vector<T>& u) {
  int n = u.n;
  int iamax;
  T amax =T(0);
  for (int i = 0; i < n; i++) {
    if (u.data[i] > amax) {
      amax = u.data[i];
      iamax = i;
    }
  }
  return iamax;
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
    u.data[i] *= a;
}

// destructive clipping
template<typename T>
inline void sc_clip(const T& a, sc_vector<T>& u) { 
  register int n = u.n;
  for (int i = 0 ; i < n; i++)  {
    T ui = u[i];
    u[i] = (ui > a)? a : ( (ui < -a)? -a : ui );
  }
}

// destructive in-place clipping
template<typename T>
inline void sc_clip(const T& a, sc_matrix<T>& A) { 
  register int n = A.n * A.m;
  for (int i = 0 ; i < n; i++)  {
    T Ai = A[i];
    A[i] = (Ai > a)? a : ( (Ai < -a)? -a : Ai );
  }
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

// copy
template<typename T>
inline void sc_copy(const sc_sparse_vector<T>& u, sc_sparse_vector<T>& v) { 
  my_assert(u.n == v.n);
  register int i = 0;
  register int n = u.nzmax;
  my_assert(u.nzmax <= v.nzmax);
  for ( ; i < n; i++)  {
    v.indexes[i] = u.indexes[i];
    v.values[i] = u.values[i];
  }
}

template<typename T>
inline void sc_scale(const T& a, sc_sparse_vector<T>& u) { 
  register int i = 0;
  register int n = u.nnz;
  for ( ; i < n; i++)  {
    u.values[i] *= a;
  }
}

//
// matrix
//
template<typename T>
inline void sc_fill(const T& a, sc_matrix<T>& A) { 
  register int i;
  register int n = A.n * A.m;
  for (i=0 ; i < n; i++)  {
    A[i] = a;
  }
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


// copy
template<typename T>
inline void sc_copy(const sc_sparse_matrix<T>& A, sc_sparse_matrix<T>& B) { 
  register int i = 0;
  my_assert(A.n == B.n);
  my_assert(A.m == B.m);
  my_assert(A.nzmax <= B.nzmax);
  int n = B.nzmax;
  // erase B first if it cannot be filled totally by A
  if (A.nzmax < B.nzmax) {
    for ( ; i < n; i++) {
      B.row_indexes[i] = -1;
    }
  }
  n = A.nzmax;
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
  register int i,j;
  register int n = A.n;
  for (j=0 ; j < n; j++)  {
    sc_sparse_vector<T> Aj;
    sc_scale(a,Aj);
  }
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
int sc_full_to_sparse(const sc_vector<T>& vf, sc_sparse_vector<T>& vs, const T& tol = T(1e-8)) 
{
  my_assert(vf.n==vs.n);
  register int i, is = 0;
  register int n=vf.n;
  for (i = 0; i < n; i++) {
    if (fabs(vf[i]) > tol) {
      if (is >= vs.nzmax) {
	std::cerr << "Error converting to sparse: too many nonzero elements (max=" << vs.nzmax << ")" << std::endl;
	break;
      }
      vs.values[is] = vf[i];
      vs.indexes[is++] = i;
    }
  }
  if (is < vs.nzmax)
    vs.indexes[is] = -1; 
  return is;
}

template<typename T> 
void sc_sparse_to_full(const sc_sparse_vector<T>& vs, sc_vector<T>& vf) 
{
  sc_fill(T(0),vf);
  for (register int i = 0; (vs.indexes[i] < 0) && (i < vs.nzmax); i++) {
    vf[vs.indexes[i]] = vs.values[i];
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

// destructive
template<typename T> 
void sc_scaled_add(const T& a, const sc_vector<T>& v, sc_vector<T>& w)
{
  my_assert(v.n == w.n);
  int n = v.n;
  for (int i = 0; i < n; i++)
    w.data[i] += a*v.data[i];
}

// destructive
template<typename T> void sc_scaled_add(const T& a, const sc_sparse_vector<T>& v, sc_vector<T>& w) {
  my_assert(v.n == w.n);
  for (int i = 0; (v.indexes[i] >= 0) && (i < v.nzmax); i++) {
    w[v.indexes[i]] += a*v.values[i];
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

template<typename T>  T sc_norm2(const sc_vector<T>& u)
{
  register T a = T(0);
  register int n = u.n;
  for (int i = 0; i < n; i++) {
    a += u[i]*u[i];
  }
  return a;
}

template<typename T>  
inline T sc_norm(const sc_vector<T>& u)
{
  return sqrt(sc_norm2(u));
}

template<typename T>  T sc_norm2(const sc_sparse_vector<T>& u)
{
  register T a = T();
  register int nzmax = u.nzmax;
  for (int i=0; (i <nzmax) && (u.indexes[i]>=0); i++) {
    a += u.values[i]*u.values[i];
  }
  return a;
}

template<typename T>  
inline T sc_norm(const sc_sparse_vector<T>& u) {
  return sqrt(sc_norm2(u));
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
  for ( ; (i < n) && (v.indexes[i] >= 0); i++) {
    a += v.values[i] * u[v.indexes[i]];
  }
  return a;
}
//
//======================================================================
//
template<typename T> 
inline T sc_dot(const sc_sparse_vector<T>& u, const sc_vector<T>& v)
{
  return sc_dot(v,u);
}
//
//======================================================================
//
template<typename T> 
inline T sc_dot(const sc_sparse_vector<T>& u, const sc_sparse_vector<T>& v)
{
  my_assert(u.n == v.n);
  register int i, j, iu = 0, iv = 0;
  int head_idx = 0; // which of the two index sets has the lead: 0->u, 1->v  
  int umax = u.nzmax;
  int vmax = v.nzmax;
  register T a(0);
  while ( (iu < umax) && (iv < vmax)  && (u.indexes[iu] >= 0) && (v.indexes[iv] >= 0) ) {
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
//
//======================================================================
//
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
//
//======================================================================
//
// dot product of a vector with the i-th column of a matrix
template<typename T>
inline T sc_dot(const sc_vector<T>& v, const sc_matrix<T>& A, int icol) { 
  return sc_dot(v,A.column(icol));
}
//
//======================================================================
//
template<typename T> 
inline T sc_dot(const sc_sparse_vector<T>& v, const sc_matrix<T>& A, int icol)
{
  return sc_dot(v,A.column(icol));
}
//
//======================================================================
//
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
//======================================================================
//
//
// matrix-vector products
//
// matrix-vector product 
template<typename T>
inline void sc_mul(const sc_matrix<T>& A, const sc_vector<T>& u, sc_vector<T>& v) { 
  sc_fill(T(0),v);
  sc_scaled_mac(1.0,A,u,v);
}
//
//======================================================================
//
// matrix-vector product with accumulation (Multiply and ACcumulate)
template<typename T>
inline void sc_scaled_mac(const T& alpha, 
			  const sc_matrix<T>& A, 
			  const sc_vector<T>& u, 
			  sc_vector<T>& v) { 
  my_assert(A.n == u.n);
  my_assert(A.m == v.n);
  register int m = A.m;
  register int n = A.n;
  register int i,j,iA;
  register T a;
  iA = 0;
  for (j=0; j < n; j++) {
    a = alpha*u[j];
    for (i=0; i < m; i++) {
      v[i] += a*A[iA++];
    }
  }
}
//
//======================================================================
//
template<typename T>
inline void sc_mul(const sc_matrix<T>& A, bool At, 
		   const sc_vector<T>& u, 
		   sc_vector<T>& v) { 
  if (~At) sc_mul(A,u,v);
  // A is considered transposed
  assert(A.m == u.n);
  assert(A.n == v.n);
  register int m = A.n;
  register int n = A.m;
  register int i,j,iA;
  register T a;
  iA = 0;
  sc_fill(T(0),v);
  for (i=0; i < m; i++) {
    v[i] = 0.0;
    for (j=0; j < n; j++) {
      v[i] += u[j]*A[iA++];
    }
  }
}
//
//======================================================================
//
template<typename T>
inline void sc_scaled_mac(const T& alpha, 
			  const sc_matrix<T>& A, 
			  const sc_sparse_vector<T>& u, 
			  sc_vector<T>& v) { 
  assert(A.n == u.n);
  assert(A.m == v.n);
  register int m = A.m;
  //register int n = A.n;
  register int i,j,iA;
  register T a;
  register int nzi;
  for (nzi=0; (nzi < u.nzmax) && (u.indexes[nzi]>=0); nzi++) {
    j = u.indexes[nzi];
    a = alpha*u.values[nzi];
    iA = m*j;
    for (i=0; i < m; i++) {
      v[i] += a*A[iA+i];
    }
  }
}

template<typename T>
inline void sc_mul(const sc_matrix<T>& A, const sc_sparse_vector<T>& u, sc_vector<T>& v) { 
  sc_fill(T(0),v);
  sc_scaled_mac(1.0,A,u,v);
}

//
//======================================================================
//
// vector-matrix product 
template<typename T>
inline void sc_mul(const sc_vector<T>& u, const sc_matrix<T>& A, sc_vector<T>& v) { 
  my_assert(A.m == u.n);
  my_assert(A.n == v.n);
  T a;
  int m = A.m;
  int n = A.n;
  int iA = 0; // automatically jump from one column to the next by adding one
  for (int j=0; j < n; j++) {
    a = T(0);
    for (int i=0; i < m; i++, iA++) {
      a += A.data[iA]*u.data[i];
    }
    v.data[j]=a;
  }
}
//
//======================================================================
//
// S  x D 
// St x D -> many times s x D
// S  x Dt
// St x Dt -> (D x S)t -> many times D x s
// D  x S -> many times D x s
// Dt x S -> St x D -> many times s x D
// D x St 
// Dt x St  -> S x D
//
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
    for (nzi=0; (nzi < u.nzmax) && ((i=u.indexes[nzi])>=0); nzi++) {
      a += A[iA+i] * u.values[nzi];
    }
    v[j]=a;
    iA += m;
  }
}
//
//======================================================================
//
//
// matrix-matrix product 
//
template<typename T>
void sc_scaled_mac(const T& alpha,
		   const sc_matrix<T>& A, bool At, 
		   const sc_matrix<T>& B, bool Bt, 
		   const T& beta,
		   sc_matrix<T>& C) { 
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
      C[iC] = beta*C[iC] + alpha*a;
      iC++;
      iA0 += dAi;
    }
    iB0 += dBj;
  }
}
//
//======================================================================
//
template<typename T>
void sc_scaled_mac(const T& alpha,
		   const sc_matrix<T>& A, bool At, 
		   const sc_sparse_matrix<T>& B, bool Bt, 
		   const T& beta,
		   sc_matrix<T>& C) {
  if (~At && Bt) { // A*B^T
    my_assert(A.n == B.n);
    my_assert(A.m == C.m);
    my_assert(B.m == C.n);
    int K = B.n;
    if (beta != 1.0) sc_scale(beta,C);
    // product computed as \sum_{k=0}^{K-1}{A_k (B_k)^T}
    for (int k=0; k < K; k++) {
      const sc_sparse_vector<T> Bkt = B.column(k);
      const sc_vector<T> Ak = A.column(k);
      sc_rank_one_update(alpha,Ak,Bkt,C);
    }
  } else if (~At && ~Bt) { // A*B
    // 
    // implemented using A*v for each column of C
    //
    my_assert(A.n == B.m);
    my_assert(A.m == C.m);
    my_assert(B.n == C.n);
    int n = C.n;
    for (int j=0; j < n;j++) {
      sc_sparse_vector<T> Bj = B.column(j);
      sc_vector<T> Cj = C.column(j);
      if (beta != 1.0) 
	    sc_scale(beta,Cj);
      sc_scaled_mac(alpha,A,Bj,Cj);
    }
  } else {
    std::cerr << "operation not implemented yet!" << std::endl; 
    exit(1);
  }
}
//
//======================================================================
//
template<typename T>
void sc_scaled_mac(const T& alpha,
		   const sc_sparse_matrix<T>& A, bool At, 
		   const sc_sparse_matrix<T>& B, bool Bt, 
		   const T& beta,
		   sc_matrix<T>& C) { 
  if (~At && Bt) { // A*B^T
    my_assert(A.n == B.n);
    my_assert(A.m == C.m);
    my_assert(B.m == C.n);
    int K = B.n;
    if (beta != 1.0) sc_scale(beta,C);
    // product computed as \sum_{k=0}^{K-1}{A_k (B_k)^T}
    for (int k=0; k < K; k++) {
      sc_sparse_vector<T> Bkt = B.column(k);
      sc_sparse_vector<T> Ak = A.column(k);
      int nnz = Bkt.nzmax;
      for (int js=0; js < nnz; js++) {
	int j = Bkt.indexes[js];
	if (j < 0) 
	  break;
	sc_vector<T> Cj = C.column(j);
	sc_scaled_add(alpha*Bkt.values[js],Ak,Cj);
      }
    }
  } else { // not implemented yet
    std::cerr << "operation not implemented yet!" << std::endl; 
  }
}
//
//======================================================================
//
template<typename T>
inline void sc_mul(const sc_matrix<T>& A, bool At, 
		   const sc_matrix<T>& B, bool Bt, 
		   sc_matrix<T>& C) {
  sc_fill(T(0),C);
  sc_scaled_mac(T(1.0),A,At,B,Bt,1.0,C);
}
//
//======================================================================
//
template<typename T>
inline void sc_mac(const sc_matrix<T>& A, bool At, 
		   const sc_matrix<T>& B, bool Bt, 
		   sc_matrix<T>& C) {
  sc_scaled_mac(T(1.0),A,At,B,Bt,1.0,C);
}
//
//======================================================================
//
// gram matrix product: G=At*A
template<typename T>
inline void sc_gram(const sc_matrix<T>& A, sc_matrix<T>& G) { 
  sc_mul(A,true,A,false,G);
}
//
//======================================================================
//
// convariance-like matrix product : S=A*At
template<typename T>
inline void sc_cov(const sc_matrix<T>& A, sc_matrix<T>& S) { 
  sc_mul(A,false,A,true,S);
}
//
//======================================================================
//
template<typename T> 
void sc_dump(const sc_matrix<T>& A, std::ostream& out) {
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
//
//======================================================================
//
template<typename T> 
void sc_print(const sc_matrix<T>& A, const char* name) {
   std::cerr << name << std::endl;
   sc_dump(A,std::cerr);
}
//
//======================================================================
//
template<typename T> 
void sc_dump(const sc_vector<T>& v, std::ostream& out)  {
   out.precision(10);
   int _n = v.n;
   for (int i = 0; i<_n; ++i) {
      out << v[i] << ' ';
   }
   out << std::endl;
};
//
//======================================================================
//
template<typename T> 
void sc_print(const sc_vector<T>& v, const char* name)  {
   std::cerr << name << std::endl;
   sc_dump(v,std::cerr);
};
//
//======================================================================
//
template<typename T> 
void sc_print(const sc_sparse_vector<T>& v, const char* name)  {
   std::cerr << name << " (sparse)" << std::endl;
   sc_dump(v,std::cerr);
};
//
//======================================================================
//
template<typename T> 
void sc_dump(const sc_sparse_vector<T>& v, 
	      std::ostream& out)  {
   out.precision(5);
   out << std::fixed;
   int i,nzi,j;
   for (i = 0,nzi = 0; i < v.n; ++i) {
     if (nzi >= v.nzmax) {
       out << T(0) << " ";
       continue;
     }
     j = v.indexes[nzi];
     if (j < 0) {
       out << T(0) << " ";
       nzi = v.nzmax + 1; // from now on all are zeros
       continue;
     } else {
       if (i == j) {
	 out << v.values[nzi] << " ";
	 nzi++;
       } else {
	 out << T(0) << " ";
       }
     }
   }
   out << std::endl;
}
//
//======================================================================
//
template <typename T> 
void sc_dump(const sc_sparse_matrix<T>& A, std::ostream& out) {
   out.precision(5);
   out << std::fixed;
   for (int i = 0; i< A.m; ++i) {
     for (int j = 0; j< A.n; ++j) {
       out << A(i,j) << " ";
     }
     out << std::endl;
   }
};

template<typename T>
void sc_print(const sc_sparse_matrix<T>& A, 
	      const char* name) {
   std::cerr << name << std::endl;
   std::cerr << A.m << " x " << A.n << " (sparse)" << std::endl;
   sc_dump(A,std::cerr);
}
#undef my_assert
#endif
