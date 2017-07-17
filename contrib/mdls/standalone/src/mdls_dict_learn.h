#ifndef MDLS_DICT_LEARN
#define MDLS_DICT_LEARN
#include "mdls_linalg.h"
#include "huber.h"
#include "mdls_pgm.h"
#include "mdls_dict_display.h"
//#include "mdls_ist.h" // for debugging only

enum { DEBUG_NONE=0, DEBUG_MINIMAL=1, DEBUG_MILD=2, DEBUG_MEDIUM=3,
       DEBUG_HEAVY=4, DEBUG_INSANE=5 };

char default_cache_dir[] = "/workspace/sparse/cache";
/**
 * stores data used for learning a dictionary
 * will be multi-class in the future. For now very simple.
 */
template<typename T>
struct csc_learning_data {
  const sc_matrix<T>& train;
  const sc_matrix<T>& test;

  csc_learning_data(const sc_matrix<T>& _train, 
		    const sc_matrix<T>& _test): 
  train(_train), test(_test) {}
};

/**
 * stores the current state of the dictionary learning algorithm,
 * including current dictionary, statistics, etc.
 */
template<typename T>
struct csc_learning_state {
  sc_matrix<T> D; // current dictionary
  double N; // efective number of samples processed
  int r; // current iteration (for resuming)

  // statistics
  sc_matrix<T> AAt;
  sc_matrix<T> XAt;
  sc_vector<T> rho; // occurence of each atom
  sc_matrix<T> xrho; // co-occurence matrix of atoms


  // cost function related  
  sc_matrix<T> maxMu; // unused for now
  sc_matrix<T> meanMu; // unused for now
  sc_vector<T> L; // codelength at each validation iteration (xval)

  csc_learning_state():
  D(), N(0), r(0), AAt(), XAt(), rho(), xrho(), maxMu(), meanMu() {}

  csc_learning_state(int M, int K, int J) {    
    XAt.m = M; XAt.n = K; XAt.data = new T[M*K];
    AAt.m = K; AAt.n = K; AAt.data = new T[K*K];

    D.m   = M; D.n   = K; D.data   = new T[M*K];

    rho.n = K;            rho.data = new T[K];
    xrho.m= K; xrho.n= K; xrho.data= new T[K*K];

    maxMu.m  = K;  maxMu.n = J; maxMu.data  = new T[K*J];
    meanMu.m = K; meanMu.n = J; meanMu.data = new T[K*J];
    L.n = J;            L.data = new T[J];
    N = 0;
    r = 0;
  }

  ~csc_learning_state() {
    if (AAt.data != NULL) delete[] AAt.data;
    if (XAt.data != NULL) delete[] XAt.data;
    if (meanMu.data != NULL) delete[] meanMu.data;
    if (maxMu.data != NULL) delete[] maxMu.data;
    if (rho.data != NULL) delete[] rho.data;
    if (xrho.data != NULL) delete[] xrho.data;
    if (D.data != NULL) delete[] D.data;
    if (L.data != NULL) delete[] L.data;
  }

};

//
//=========================================================================
//

/**
 * stores the dictionary learning algorithm parameters
 */
struct csc_learning_params {
  int max_iter;
  int min_change;
  int batch_size;
  int test_size;
  double inertia;
  double discard_unused_atoms; // not used yet
  double discard_const_patches; // not used yet
  int xval_step; // every how many batches do we validate
  double mu0; // incoherence term, not used yet
  char* cache_dir; // where to backup the state of the algorithm
  int debug; // debug level
  bool dump; // dump intermediate dictionaries as pgm files

  csc_learning_params():
    max_iter(100), min_change(1e-4), batch_size(0), test_size(0),
    inertia(0.8),discard_unused_atoms(0), discard_const_patches(1e-1),
      xval_step(10), mu0(0), cache_dir(default_cache_dir), debug(DEBUG_NONE),
      dump(false)
  {
  }
};

//
//=========================================================================
//

/**
 * Fast incremental dictionary learning using CSC-MP for sparse coding
 * Finds a local minimum of
 * 
 * min_{D,A}  L(X-DA) + L(A) s.t. ||d_k|| <= 1
 * 
 * where L(D,A) is the quadratic loss corresponding to a 
 * Gaussian prior on the approximation error.
 *
 * X ........... training data, given as dense MxN matrix
 * params ...... learning algorithm parameters
 * state ....... on entry, initial state. On exit, final state.
 */
template<typename T> 
void mdls_dict_learn_csc_l2(const csc_learning_data<T>& data, // learn input
			    const csc_coding_params& cparams, // input
			    const csc_learning_params& mparams, // input
			    csc_learning_state<T>& state); // input/output

//
//=========================================================================
//
/**
 * Fast incremental dictionary learning using CSC-MP for sparse coding
 * Finds a local minimum of
 * 
 * min_{D,A}  L(X-DA) + L(A) s.t. ||d_k|| <= 1
 * 
 * where L(D,A) is the Huber loss, which is more consistent than L2 for the
 * observed empirical statistics of approximation errors.
 *
 * X ........... training data, given as dense MxN matrix
 * params ...... learning algorithm parameters
 * state ....... on entry, initial state. On exit, final state.
 */
template<typename T> 
void mdls_dict_learn_csc_huber(const csc_learning_data<T>& data, // learn input
				const csc_coding_params& cparams, // input
				const csc_learning_params& mparams, // input
				csc_learning_state<T>& state); // input/output


//
//=========================================================================
//
/**
 * Minimal efficient dictionary update. 
 * Finds the minimum of 
 * 
 * min_D  ||X-DA||_2^2 s.t. ||d_k|| <= 1
 *
 * using block coordinate descent. This is the algorithm
 * described in Mairal at NIPS 2009. 
 *
 * As this is used in an alternating descent problem, finding the actual
 * minimum is usually not necessary and we may be happy with a single
 * iteration.
 *
 * AAt ......... A*A^T statistic (KxK)
 * XAt ......... X*A^T statistic (MxK)
 * J ........... number of iterations (1)
 * D ........... dictionary to be updated
 * aux ......... auxiliary vector for internal computations (1xK)
 */
template<typename T> 
void mdls_dict_update_l2(const sc_matrix<T>& AAt, 
			 const sc_matrix<T>& XAt, 
			 const int J,
			 sc_matrix<T>& D, 
			 sc_vector<T>& aux);

//
//=========================================================================
//
/**
 * Dictionary update for hyperbolic (quasi-l1) loss
 * Finds the minimum of 
 * 
 * min_D  L(X-DA) s.t. ||d_k|| <= 1
 *
 * where L(x) = sqrt(x^2+e^2), where e is a small constant.
 *
 * using block coordinate descent. This is the algorithm
 * described in Mairal at NIPS 2009. 
 *
 * AAt ......... A*A^T statistic (KxK)
 * XAt ......... X*A^T statistic (MxK)
 * J ........... number of iterations (1)
 * D ........... dictionary to be updated
 * aux ......... auxiliary vector for internal computations (1xK)
 */
template<typename T> 
void mdls_dict_update_hyperbolic(const sc_matrix<T>& AAt, 
		      const sc_matrix<T>& XAt, 
		      const int J,
		      sc_matrix<T>& D, 
				 sc_vector<T>& aux);


//
//=========================================================================
//
/**
 * Dictionary update for Huber loss (approximates LG distribution)
 * Finds the minimum of 
 * 
 * min_D  L(X-DA) s.t. ||d_k|| <= 1
 *
 * where L(x) = Huber(x) = (|x| < d)? x^/2 : d(|x|-d/2)
 *
 * The Hessian is estimated
 * using the BFGS method (starting with B=I). 
 * Unfortunately this function cannot be made to depend only
 * on the summary statistics AAt and XAt, so it cannot
 * be implemented online (maybe yes, if we could argue
 * that as the dictionary converges, the past residuals will
 * not change, some stochastic descent argument I guess...
 *
 * R ........... current residual: R=X-DA (DESTROYED)
 * A ........... coefficients
 * d ........... Huber funcion parameter
 * J ........... number of iterations (1)
 * D ........... dictionary to be updated
 * Hinv ........ storage for inverse of Hessian (KxK)
 * aux ......... auxiliary vector for internal computations (1xK)
 */
template<typename T> 
void mdls_dict_update_huber(const sc_matrix<T>& X, 
			    const sc_sparse_matrix<T>& A, 
			    const T d,
			    const int J,
			    sc_matrix<T>& D, 
			    sc_huber_workspace<T>& wks);

//
//=========================================================================
//
/**
 * Dictionary update for generic loss specified by lookup table.
 * 
 * min_D  L(X-DA) s.t. ||d_k|| <= 1
 *
 * where L(x) is given as a lookup table {x_i,l_i}. This
 * function finds L'(x) using the GSL cubic interpolation functions,
 * for a resulting gradient L'(X-DA) * A^T, and estimates the Hessian
 * using the BFGS method (starting with B=I). See parametric_lut.h
 *
 * AAt ......... A*A^T statistic (KxK)
 * XAt ......... X*A^T statistic (MxK)
 * J ........... number of iterations (1)
 * D ........... dictionary to be updated
 * aux ......... auxiliary vector for internal computations (1xK)
 */
template<typename T> 
void mdls_dict_update_lookup(const sc_matrix<T>& AAt, 
			     const sc_matrix<T>& XAt,
			     const parametric_lut& lut,
			     const int J,
			     sc_matrix<T>& D, 
			     sc_vector<T>& aux);
//
//=========================================================================
//=========================================================================
//=========================================================================
//
#include <ctime>

template<typename T> 
void mdls_random_subsample(const sc_matrix<T>& sample,
			   sc_matrix<T>& subsample) {
  assert(sample.m == subsample.m);
  for (int destidx = 0; destidx < subsample.n; destidx++) {
    int srcidx = (int) (sample.n*double(rand())/double(RAND_MAX));    
    sc_vector<double> srcX = sample.column(srcidx);
    sc_vector<double> destX  = subsample.column(destidx);
    sc_copy(srcX, destX);
  }
}

//
//=========================================================================
//=========================================================================
//=========================================================================
//
template<typename T> 
double mdls_avg_sparsity(const sc_sparse_matrix<T>& A) {
  double nnz = T(0);
  for (int j=0; j < A.n; j++) {
    sc_sparse_vector<T> Aj = A.column(j);
    nnz += (T) Aj.nnz();
  }
  return nnz/double(A.n);
}

//
//=========================================================================
//=========================================================================
//=========================================================================
//
template<typename T> 
void mdls_dict_learn_csc_l2(const csc_learning_data<T>& data, // learn input
			    const csc_coding_params& cparams, // input
			    const csc_learning_params& mparams, // input
			    csc_learning_state<T>& state) // input/output
{
  int M  = data.train.m;
  int N = data.train.n;
  int Nb = mparams.batch_size;
  int Nt = data.test.n;
  int nzmax = cparams.nzmax;
  int K  = state.D.n; 
  int r = state.r;
  int J = mparams.max_iter;
  double thres = mparams.min_change;
  time_t t0,t1;
  time(&t0);
  //
  // allocate working space
  // 
  sc_sparse_matrix<T> Ab; // stores sparse code for current batch
  sc_matrix<T> Xb; // stores subsample for current batch
  sc_vector<T> Lb;  // store codelength obtained for each sample in batch

  sc_sparse_matrix<T> At; // stores sparse code for testing samples
  sc_vector<T> Lt;  // store codelength obtained for  testing samples
  sc_vector<T> aux;

  if (Nb > 0) { // for mini-batch mode, allocate space for sub sample
    Xb.allocate(M,Nb);
    Ab.allocate(K,Nb,nzmax);
    Lb.allocate(Nb);
  } else { 
    Xb = data.train;
    Ab.allocate(K,N,nzmax);
    Lb.allocate(N);
  }
  At.allocate(K,Nt,nzmax);
  Lt.allocate(Nt);

  aux.allocate(K); // auxiliary for dictionary update
  //
  // for IST (debug)
  //
  //sc_vector<T> W;
  //W.allocate(K);
  //sc_fill(255.0,W);
  //ist_workspace<T> ist_wks;
  //ist_wks.allocate(K,nzmax);
  //ist_params iparams;
  //iparams.c = 30.0;
  //iparams.maxIter = 2000;
  //
  // initial cost
  //
  int stuck = 0;
  csc_fss_lg_l(data.test, state.D, cparams, At, Lt);
  //mdls_ist(data.test, state.D, W, iparams, At, ist_wks); // REMOVE ME
  double cost = sc_sum(Lt)/((double)Nt);
  double annz = mdls_avg_sparsity(At);
  std::cout << "L=" << cost << " <nnz>=" << annz << std::endl;
  
  //
  // main loop
  //  
  char dicfname[64];
  FILE* dicf = NULL;
  sc_matrix<double> I;
  if (mparams.dump) {
    mdls_dict_display(state.D,I);
    snprintf(dicfname,63,"D_l2_i%03d.pgm",0);
    dicf = fopen(dicfname,"w");
    mdls_pgm_write(I.data,I.m,I.n,dicf);
    fclose(dicf);
  }
  sc_fill(T(0),state.XAt);
  sc_fill(T(0),state.AAt);
  while(r < J) {
	r++;
    //
    // 1. subsample learning data at random
    //
	if (mparams.debug >= DEBUG_MEDIUM)
		std::cout << "Subsampling." << std::endl;
    if (Nb > 0) { 
      mdls_random_subsample(data.train, Xb);
    }
    //
    // 2. perform sparse coding 
    //
	if (mparams.debug >= DEBUG_MEDIUM)
		std::cout << "Coding." << std::endl;
    csc_fss_lg_l(Xb, state.D, cparams, Ab, Lb);
    // HACK: ver que pasa con IST
    //mdls_ist(Xb, state.D, W, iparams, Ab, ist_wks);

    //
    // 3. update statistics
    // 
    // XAt += Xb*(Ab)^T
	if (mparams.debug >= DEBUG_MEDIUM)
		std::cout << "Updating statistics." << std::endl;
    double oldN = state.N*mparams.inertia;
    double newN = oldN + ((Nb >0)? Nb : N);
    sc_scaled_mac(1.0/newN,Xb,false,Ab,true,oldN/newN,state.XAt);
    sc_scaled_mac(1.0/newN,Ab,false,Ab,true,oldN/newN,state.AAt);
    state.N = newN;
    //sc_scaled_mac(1.0,Xb,false,Ab,true,1.0,state.XAt);
    //sc_scaled_mac(1.0,Ab,false,Ab,true,1.0,state.AAt);
    //
    // 4. update dictionary
    //
	if (mparams.debug >= DEBUG_MEDIUM)
		std::cout << "Updating dictionary." << std::endl;
    mdls_dict_update_l2(state.AAt,
			state.XAt,
			1, // one iteration
			state.D,
			aux);
    //
    // 5. evaluate progress/stopping criterion 
    //    only if cross-validation is performed in this iteration 
    //
    if ((r % mparams.xval_step) == 0) { // is this is an xval iteration?
    	if (mparams.debug >= DEBUG_MEDIUM)
    		std::cout << "Cross validating." << std::endl;
      //
      // optionally dump current dictionary to file
      //
      if (mparams.dump) {
	mdls_dict_display(state.D,I);
	snprintf(dicfname,63,"D_l2_i%03d.pgm",r);
	dicf = fopen(dicfname,"w");
	mdls_pgm_write(I.data,I.m,I.n,dicf);
	fclose(dicf);
      }
      //
      // 5a. sparse coding of testing set
      //
      csc_fss_lg_l(data.test, state.D, cparams, At, Lt);
      double new_cost = sc_sum(Lt)/double(Nt);
      double annz = mdls_avg_sparsity(At);
      std::cout << "L=" << new_cost << " <nnz>=" << annz << std::endl;
      //
      // 5b. average codelength obtained
      //
      //
      // 5c. did it improve significantly?
      // (signed comparison means we also stop if codelength increased)
      //
      if ((cost - new_cost) < thres) {
	stuck++;
	if (stuck >= 3) {
	  break; // no further improvement in last 3 iterations
	}
      } else {
	stuck = 0;
      }
      cost = new_cost;
    }
    time(&t1);    
  } // main loop over batches
	if (mparams.debug >= DEBUG_MEDIUM)
		std::cout << "Finishing." << std::endl;
  
  //
  // free working space
  //W.free();
  //ist_wks.free();
  //
  if (Nb > 0) 
    Xb.free();
  Ab.free();
  Lb.free();

  At.free();
  Lt.free();
  aux.free();
}

//
//=========================================================================
//=========================================================================
//=========================================================================
//
//
//=========================================================================
//=========================================================================
//=========================================================================
//

template<typename T> 
void mdls_dict_learn_csc_huber(const csc_learning_data<T>& data, // learn input
			       const csc_coding_params& cparams, // input
			       const csc_learning_params& mparams, // input
			       csc_learning_state<T>& state) // input/output
{
  int M  = data.train.m;
  int N = data.train.n;
  int Nb = mparams.batch_size;
  int Nt = data.test.n;
  int nzmax = cparams.nzmax;
  int K  = state.D.n; 
  int r = state.r;
  int J = mparams.max_iter;
  double thres = mparams.min_change;
  time_t t0,t1;
  time(&t0);
  //
  // allocate working space
  // 
  sc_sparse_matrix<T> Ab; // stores sparse code for current batch
  sc_matrix<T> Xb; // stores subsample for current batch
  sc_vector<T> Lb;  // store codelength obtained for each sample in batch

  sc_sparse_matrix<T> At; // stores sparse code for testing samples
  sc_vector<T> Lt;  // store codelength obtained for  testing samples
  sc_vector<T> aux;

  if (Nb > 0) { // for mini-batch mode, allocate space for sub sample
    Xb.allocate(M,Nb);
    Ab.allocate(K,Nb,nzmax);
    Lb.allocate(Nb);
  } else { 
    Xb = data.train;
    Ab.allocate(K,N,nzmax);
    Lb.allocate(N);
  }
  At.allocate(K,Nt,nzmax);
  Lt.allocate(Nt);

  aux.allocate(K); // auxiliary for dictionary update
  //
  // initial cost
  //
  csc_fss_lg_l(data.test, state.D, cparams, At, Lt);
  double cost = sc_sum(Lt)/double(Nt);
  double annz = mdls_avg_sparsity(At);
  std::cout << "L=" << cost << " <nnz>=" << annz << std::endl;

  char dicfname[64];
  sc_matrix<double> I;
  FILE* dicf;
  if (mparams.dump) {
    mdls_dict_display(state.D,I);
    snprintf(dicfname,63,"D_huber_i%03d.pgm",0);
    dicf = fopen(dicfname,"w");
    mdls_pgm_write(I.data,I.m,I.n,dicf);
    fclose(dicf);
  }
  //
  // main loop
  //
  int stuck = 0;
  sc_huber_workspace<T> huber_wks;
  huber_wks.allocate(M);
  //
  // huber parameter
  //
  //T h = T(5.0)*cparams.sigma;
  T h = T(5.0)*cparams.sigma;
  sc_fill(T(0),state.XAt);
  sc_fill(T(0),state.AAt);

  while( r < J ) {
	r++;
    //
    // 1. subsample learning data at random
    //
    if (Nb > 0) { 
      mdls_random_subsample(data.train, Xb);
    }
    //
    // 2. perform sparse coding 
    //
    csc_fss_lg_l(Xb, state.D, cparams, Ab, Lb);
    //
    // 3. update statistics
    // 
    // XAt += Xb*(Ab)^T
    double oldN = state.N*mparams.inertia;
    double newN = oldN + ((Nb >0)? Nb : N);
#define V2 1
#ifdef V2
    sc_scaled_mac(1.0,state.D,false,Ab,false,-1.0,Xb);
    // here Xb = DAb-Xb
    sc_clip(h,Xb);
#endif
    // here Xb = clip(Xb-DAb)
    sc_scaled_mac(1.0/newN,Xb,false,Ab,true,oldN/newN,state.XAt);
    // now XAt stores clip(X-DA)*A^T
    sc_scaled_mac(1.0/newN,Ab,false,Ab,true,oldN/newN,state.AAt);
    state.N = newN;
    // here Xb = clip(Xb-DAb)
    //sc_scaled_mac(1.0,Xb,false,Ab,true,1.0,state.XAt);
    // now XAt stores clip(X-DA)*A^T
    //sc_scaled_mac(1.0,Ab,false,Ab,true,1.0,state.AAt);
    //
    // S(t) = (N*S(-1)+s)/(N+n) = N/(N+n)S(-1) + 1/n
    // 4. update dictionary
    //
    int du_iter = 1;
    //T huber_par = T(1.0)*cparams.sigma;
    // stochastic-like descent: 
    // update is only based on Xb and Ab
    // however, BFGS hessian estimation is cumulative
    //
#ifdef V2
    mdls_dict_update_huber_v2(state.AAt,
			      state.XAt, 
			      h, 
			      du_iter,
			      state.D, 
			      huber_wks);
#else
    mdls_dict_update_huber(Xb,Ab, h, du_iter,
			   state.D, huber_wks);
#endif
    //
    // 5. evaluate progress/stopping criterion 
    //    only if cross-validation is performed in this iteration 
    //
    if ((r % mparams.xval_step) == 0) { // is this is an xval iteration?
      if (mparams.dump){ 
	mdls_dict_display(state.D,I);
	snprintf(dicfname,63,"D_huber_i%03d.pgm",r);
	dicf = fopen(dicfname,"w");
	mdls_pgm_write(I.data,I.m,I.n,dicf);
	fclose(dicf);
      }
      //
      // 5a. sparse coding of testing set
      //
      csc_fss_lg_l(data.test, state.D, cparams, At, Lt);
      double new_cost = sc_sum(Lt)/double(Nt);
      double annz = mdls_avg_sparsity(At);
      std::cout << "L=" << new_cost << " <nnz>=" << annz << std::endl;
      //
      // 5c. did it improve significantly?
      // (signed comparison means we also stop if codelength increased)
      //
      if ((cost - new_cost) < thres) {
	stuck++;
	if (stuck >= 3) {
	  break; // no further improvement in last 3 iterations
	}
      } else {
	stuck = 0;
      }
      cost = new_cost;
    }
    time(&t1);    
  } // main loop over batches
  //
  // free working space
  //
  huber_wks.free();
  if (Nb > 0) 
    Xb.free();
  Ab.free();
  Lb.free();

  At.free();
  Lt.free();
  aux.free();
}

//
//=========================================================================
//=========================================================================
//=========================================================================
//

template<typename T> 
void mdls_dict_update_l2(const sc_matrix<T>& AAt, 
		      const sc_matrix<T>& XAt, 
		      const int J,
		      sc_matrix<T>& D, 
		      sc_vector<T>& aux)
{ 
  int K = D.n;
  for (int i = 0; i < J; i++) { // iteration
    for (int k=0; k < K; k++) { // atoms
      T AAkk = AAt(k,k);
      double C = (AAkk > 1e-5)? 1.0/AAkk: 1e5; //damped diagonalized Newton
      //
      // 1) d_k <- d_k + (1/AAt(k,k))*( XAt(:,k)-D*AAt(:,k) )
      //
      sc_vector<T> dk = D.column(k);
      sc_vector<T> XAtk = XAt.column(k);
      sc_vector<T> AAtk = AAt.column(k);
      sc_scaled_add(C,XAtk,dk); // d_k += XAt(:,k)/AAt(k,k)
      sc_scaled_mac(C,D,AAtk,dk); // d_k += D*AAt(:,k)/AAt(k,k)
      //
      // 2) d_k <- (1/max{||d_k|,1})*d_k
      //
      C = sc_norm(dk);
      if (C >= 1.0) {
    	  sc_scale(1.0/C,dk);
      }
    } // for each atom
  } // iteration
} // function

template<typename T> 
void mdls_dict_update_hyperbolic(const sc_matrix<T>& AAt, 
		      const sc_matrix<T>& XAt, 
		      const int J,
		      sc_matrix<T>& D, 
				 sc_vector<T>& aux);

//
//=========================================================================
//=========================================================================
//=========================================================================
//
//
// simplified version
//
template<typename T> 
void mdls_dict_update_huber_v2(const sc_matrix<T>& AAt, 
			       const sc_matrix<T>& dD,       
			       const T d,
			       const int J,
			       sc_matrix<T>& D, 
			       sc_huber_workspace<T>& wks)
{
  int K = D.n;
  for (int i = 0; i < J; i++) { // iteration
    for (int k=0; k < K; k++) { // atoms
      T AAkk = AAt(k,k);
      assert(!isnan(AAkk));
      double C = (AAkk > 1e-5)? 1.0/AAkk: 1e5; //damped diagonalized Newton
      sc_vector<T> dk  = D.column(k);
      sc_vector<T> dDk = dD.column(k);
      sc_scaled_add(-C,dDk,dk); // Dk -< Dk + sk*Gk
      C = sc_norm(dk);
      assert(!isnan(C));
      if (C >= 1.0) {
	    sc_scale(1.0/C,dk);
      }
    } // for each atom
  } // iteration
}

template<typename T> 
void mdls_dict_update_huber(const sc_matrix<T>& X, 
			    const sc_sparse_matrix<T>& A, 
			    const T d,
			    const int J,
			    sc_matrix<T>& D, 
			    sc_huber_workspace<T>& wks)
{
  int K = D.n;
  int N = X.n;
  double thres = 1e-6*N;
  for (int iter = 0; iter < J; iter++) { // outer: number of BCD iterations
#ifdef DEBUG
    std::cout << "inner iter " << iter << std::endl;
#endif
    for (int k=0; k < K; k++) { // middle: for each atom
      sc_vector<T> Dk = D.column(k);      
      // store previous atom
      //
      // block update: projected gradient
      //
      int atom_iter = 0;
#ifdef DEBUG
      std::cout << "atom " << k << std::endl;
#endif
      while (atom_iter < 1) {
	//
	// 1. evaluate partial gradient for current atom k, 
	// Gk = [Huber'(X-DA;d)]*A^k = Clip(X-DA;d)*A^k
	//
	atom_iter++;
	sc_fill(T(0),wks.Gk);
	T AAkk = T(0);
	for (int j=0; j < N; j++) {
	  sc_sparse_vector<T> Aj = A.column(j);
	  T Akj = Aj[k];
	  if (Akj == T(0)) { // skip if atom not used by this sample
	    continue;
	  } else {
	    AAkk += Akj*Akj;
	    sc_vector<T> Xj = X.column(j);
	    sc_copy(Xj,wks.Rj);
	    sc_scaled_mac(-1.0,D,Aj,wks.Rj); // Rj = X_j - D*A_j
	    //std::cout << "Akj=" << Akj << std::endl;
	    //std::cout << "||Rj||_2=" << sc_norm(wks.Rj) << std::endl;
	    sc_clip(d,wks.Rj); // clipping function 
	    //std::cout << "||clip(Rj)||_2=" << sc_norm(wks.Rj) << std::endl;
	    sc_scaled_add(Akj,wks.Rj,wks.Gk);
	    //std::cout << "||Gk||_2=" << sc_norm(wks.Gk) << std::endl;
	  } 
	} // for each data sample
	if (AAkk == 0)
	  continue;
	//
	//	sc_print(wks.Gk,"Gk");
	//
	// 2. find update step using fixed step
	//    Armijo may be too slow here!
	//   also I don't know how to find a good scaling matrix
	//   for this case, so I use the L2 diagonal approximation
	//   as if it were an L2 update, divided by 10, to counter
	//   for the deviation in the model... really ugly hack
	//   but at least it is a valid descent direction...
	//T s(0), beta(0.25), alpha(0.1);
	T sk = 1.0/AAkk;
	sc_copy(Dk,wks.prevDk);
	sc_scaled_add(sk,wks.Gk,Dk); // Dk -< Dk - sk*Gk
	T nk = sc_norm(Dk);
	if (nk > 1.0) 
	  sc_scale(1.0/nk,Dk); // normalize Dk
	//
	// check for convergence
	//
	T dD = sc_dist(Dk,wks.prevDk);
#ifdef DEBUG
	if ((atom_iter % 10) == 0 )
	  std::cout << "||dD||_2=" << dD << std::endl;
#endif
	if (dD < thres) {
#ifdef DEBUG
	  std::cout << "Converged." << std::endl;
#endif
	  break;
	}
      } // BFGS iteration for a single atom
    } //  for each atom
  } // for each outer loop of the block-coordinate descent
 } // function


//=========================================================================
//=========================================================================
//=========================================================================

#endif
