#ifndef MDLS_IST_H
#define MDLS_IST_H
#include <cmath>
#include "mdls_linalg.h"
#include "mdls_shrink.h"

template<typename T>
struct ist_workspace {

  sc_matrix<T> cDtDpI; // gram matrix of dictionary
  sc_vector<T> cDtXj; // stores current correlation between D and j-th sample
  sc_vector<T> Afull; // temporary full representation of current iterate
  sc_sparse_vector<T> Acurr; // current iterate
  sc_sparse_vector<T> Aprev; // previous iterate
  sc_sparse_vector<T> Aprev2; // second prev. iterate

  ist_workspace(int K, int nzmax) {
    allocate(K,nzmax);
  }

  void allocate(int K, int L) {
    cDtDpI.allocate(K,K);
    cDtXj.allocate(K);
    Afull.allocate(K);
    Acurr.allocate(K,L);
    Aprev.allocate(K,L);
    Aprev2.allocate(K,L);
  }
  
  ist_workspace():cDtDpI(),cDtXj(),Afull(),Acurr(),Aprev(),Aprev2() {}

  ~ist_workspace() {
    free();
  }

  void free() {
    if (cDtDpI.data != NULL) 
      cDtDpI.free();
    if (cDtXj.data != NULL) 
      cDtXj.free();
    if (Afull.data != NULL) 
      Afull.free();
    if (Acurr.values != NULL)
      Acurr.free();
    if (Aprev.values != NULL) 
      Aprev.free();
    if (Aprev2.values != NULL) 
      Aprev2.free();
  }
};

/* still unused, but is a good idea */
struct ist_params {
  double alpha;
  double beta;
  double c;
  int maxIter;
  double tol;
  int skip;

  ist_params():alpha(1.0),beta(1.0),c(1.0),maxIter(200),tol(1e-6),skip(0){}
};

#define INITIAL_FULL_ITERS 3
//
// TODO
// -increasing skip size: empirically, support does not change after a few dozen iterations (40-50)
// ... or, do 'support monitoring': if support does not change for L iterations, increase skip geometrically
//
template<class T> 
void mdls_ist_precomp(const sc_matrix<T>& X, 
		      const sc_matrix<T>& D, 
		      const sc_vector<T>& W,
		      const ist_params& params,
		      sc_sparse_matrix<T>& A,  // output
		      ist_workspace<T>& wspace) // working data
{
   /* auxiliary */
   register int i,j,k;
   register T g;
   double dA;
   int M,N,K;
   double average_iters = 0;
   int base_skip = params.skip;
   int maxIter   = params.maxIter;
   T c = T(params.c);
   //T a = T(params.alpha);
   //T b = T(params.beta);
   T tol = T(params.tol);

   if (base_skip < 1)
     base_skip = 1;
   //int dsup, csup, psup, fsup;
   N = X.n;
   K = D.n;
   M = X.m;
   //sc_scale( T(1.0/(2.0*(double)c)) , W);
   
#ifdef DEBUG
          std::cout << " skip=" << base_skip << " alpha=" << alpha << " beta=" << beta << " c=" << c << std::endl;
          std::cout << "D.n=" << D.n << " D.m=" << D.m  << " cDtXj.n=" << wspace.cDtXj.n << std::endl;
#endif
   for (j = 0; j < N; j++) { // Acurr over samples
     int skip = base_skip;
     sc_vector<T> Xj = X.column(j);
     sc_sparse_vector<T> Aj = A.column(j);
     // would this be better?
     //     if (j > 1)
     //       sc_copy(A.column(j-1),Aj);

     sc_copy( Aj, wspace.Acurr );
     sc_mul( Xj, D, wspace.cDtXj ); // D^t*Xj = Xj^t*D    
     sc_scale( T(1.0/c), wspace.cDtXj );     
     for (i = 0; i < maxIter ; i++) { // IS iterations
       //dsup = 0;
        //sc_copy(wspace.Aprev,wspace.Aprev2);
       //%sc_copy(wspace.Acurr,wspace.Aprev);
       if ((i<INITIAL_FULL_ITERS) || ((i % skip) == 0)) { 
	 //
	 // full update
	 //
	 // p = (1/c)*D^t*X -(1/c)*D^t*D*A + A 
	 dA = T();
	 sc_sparse_to_full(wspace.Acurr,wspace.Afull);
	 for (k=0; k < K; k++) {
	   // preshrink : (1/c)<X,d_k> 
	   // c*Dt*D*A = At * c*(Dt*D)t = At * c*Dt*D
	   //psup = (wspace.Afull[k] != T());
	   g = wspace.cDtXj[k] - sc_dot(wspace.Acurr, wspace.cDtDpI,k);
	   T t = 0.5*W[k]/c;
	   if (g > t) {
	     g -= t;
	   } else if (g < -t) {
	     g += t;
	   } else {
	     g = T();
	   }
	   t = g - wspace.Afull[k];
	   dA += t*t;
	   //csup = (g != T());
	   wspace.Afull[k] = g;
	   //if (csup != psup) 
	   //  dsup++;
	 }	 
	 dA = sqrt(dA);
	 int nnz = sc_full_to_sparse(wspace.Afull,wspace.Acurr);
	 if (nnz >= wspace.Acurr.nzmax)
	   break;
	 if (dA <= tol) {
	   break;
 	 }
	 //
	 // next skip is twice as large
 	 //
	 if (i % skip)
	   skip = skip * 2;
       } else {
	 //
	 // partial update
	 //
	 // do computations only on current active set (indexes in currIters)
	 // element by element
	 register size_t nzi,nzn;
	 dA = T();
	 nzn = wspace.Acurr.nzmax;
	 for (nzi=0; nzi < nzn; nzi++) {
	   k = wspace.Acurr.indexes[nzi];
	   if (k-- < 1) // no more nonzero elements
	     break;
	   // preshrink : (1/c)<X,d_k> 
	   g = wspace.cDtXj[k] - sc_dot(wspace.Acurr, wspace.cDtDpI,k);
	   T t = W[k];
	   if (g > t) {
	     g -= t;
	   } else if (g < -t) {
	     g += t;
	   } else {
	     g = T();
	   }
	   wspace.Acurr.values[nzi] = g;
	 }
       }
     } // for each IST iteration 
     average_iters += i;
     sc_copy(wspace.Acurr,Aj);
   } // for each sample
   //#ifdef DEBUG
   std::cout << "average iterations =" << average_iters/(double)N << std::endl;
   //#endif   
}

/**
 * X ............ data (n x N)
 * D ............ dictionary (n x K)
 * a ............ twist constant [1.0]
 * b ............ twist constant [1.0]
 * c ............ scale factor
 * W ............ regularization parameter vector (1 x K)
 * maxIter ...... maximum number of iterations
 * tolerance .... threshold for the relative variation of the Frobenius norm
 *                of A
 * skip ......... skip update of support for this many iterations [1]
 *
 * Output:
 * A .. ......... new value for the coefficients (K x N)
 * wspace ....... IST workspace
 */
template<class T> 
void mdls_ist(const sc_matrix<T>& X, 
	      const sc_matrix<T>& D, 
	      const sc_vector<T>& W,  
	      const ist_params& params,
	      sc_sparse_matrix<T>& A,  // output
	      ist_workspace<T>& wspace) { // working data
  sc_gram(D,wspace.cDtDpI);
  mdls_ist_precomp(X, D, W, params, A, wspace);
}

#endif
