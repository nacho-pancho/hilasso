#ifndef MDLS_CSC_FSS
#define MDLS_CSC_FSS
#include "mdls_linalg_types.h"
#include "csc_types.h"
#include <iomanip>

//#define DEBUGMP

#ifdef _OPENMP // parallel 
  #include <omp.h>
#endif

/**
 * Codelength-based Sparse Coding  
 * Forward Selection (CSC-FS),
 * Laplacian+Gaussian (LG) prior on residual error and
 * Laplacian (L) prior on reconstruction coefficients.
 * Support is described using enumerative code for known N
 *
 * This is the main external interface:
 *
 * D ........ [MxK] matrix storing the dictionary
 * X ........ [MxN] matrix with data to be encoded
 * params ... csc_coding_params structure with the configuration of the algo *
 * OUTPUT
 * A ........ [KxN] sparse matrix to be filled with coefficients, has to be
able to store NZMAXxN nonzero elements (NZMAX is specifyed in params)
 * L ....... [1xN] vector to store the resulting codelength of each sample
 */
template<class T>
void csc_fss_lg_l(const sc_matrix<T>& X,
		  const sc_matrix<T>& D, 
		  const csc_coding_params& params,
		  sc_matrix<T>& A,
		  sc_vector<T>& L);

//
//=========================================================================
//
/**
 * internal structure with constant parmeters and precalculated
 * constant magnitudes of the algorithms.
 */
template<typename T>
class csc_constants {
public:
  double qx; // signal quantization step
  double qa; // coefficients quantization step  
  double maxx;   // x takes values in [0,maxx]
  double sigma; // variance of Gaussian noise in model
  bool dump_path;
  int step_mode; // 0: non-greedy, stepwise, 1: greedy matching pursuit
  int nzmax; // forced maximum number of nonzero elements in reconstruction
  int Loffset; // additional fixed codelength to be added to the final result
  parametric_lut* Lr_lut; // optional lookup table for codelength of residual r
  parametric_lut* La_lut; // optional lookup table for codelength of coeffs.  a
  sc_matrix<T> D; // copy of dictionary, used for parallel impl.
  sc_matrix<T> DtD; // Gram matrix of D, scaled by quantization step qa
  int M,K;
  int j; // current sample index, used when dumping path

  csc_constants(const sc_matrix<T>& _D, const csc_coding_params& params) 
  { 
    K = _D.n;
    M = _D.m;
    qa = params.qa;
    qx = params.qx;
    maxx = params.maxx;
    sigma = params.sigma;
    nzmax = params.nzmax;
    dump_path = params.dump_path;
    step_mode = params.step_mode;
    Loffset = params.Loffset;
    D   = sc_matrix<T>(new T[M*K],M,K);
    sc_copy(_D,D);
    DtD = sc_matrix<T>(new T[K*K],K,K);
    sc_gram(D,DtD);
  }
  
  ~csc_constants() {
    delete[] DtD.data;
    delete[] D.data;
  }
};

//
//=========================================================================
//
/**
 * easier interface using sc_minimal.h structures.
 */
template<class T>
void csc_fss_lg_l(const sc_matrix<T>& X,
		  const sc_matrix<T>& D, 
		  const csc_coding_params& params,
		  sc_sparse_matrix<T>& A,
		  sc_vector<T>& L)
{
  //
  // parallel processing preamble: find out and set number of threads
  //
  int NTHREADS = 1;
#ifdef _OPENMP
  NTHREADS =  omp_get_num_procs();
  //std::cout << "threads: " << NTHREADS << std::endl;
  omp_set_num_threads(NTHREADS);
  omp_set_nested(0);
  omp_set_dynamic(0);
#endif  
  //
  // initialization
  //  
  int M = X.m;
  int N = X.n;
  int K = D.n;
  int nzmax = params.nzmax;

  assert(X.m == D.m);

  if (A.m == 0) {
    A.allocate(K,N,nzmax);
  }
  if (L.data == NULL) {
    L.allocate(N);
  }

#if DEBUG
  //std::cout << "FSS: M=" << M << " N=" << N << " K=" << K << std::endl;
  //  std::cout << "FSS: creating lookup table..." << std::endl;
#endif

#if DEBUG
  //std::cout << "qa=" << params.qa << " qx=" << params.qx << " maxx=" << params.maxx << " sigma=" << params.sigma << std::endl;
  //std::cout << "precomputing state" << std::endl;
#endif

  //
  // initialize data structures
  //
  // constant: parameters of the algorithm (and derived constants such as DtD)
  //
  csc_constants<T>** constants = new csc_constants<T>*[NTHREADS]; //constant parameters
  csc_state<T>** state = new csc_state<T>*[NTHREADS];   //state of algorithm

  for (int t=0; t < NTHREADS; t++) {
    constants[t] = new csc_constants<T>(D,params);
    parametric_lut* lut = create_lg_parametric_lut(params.maxx, params.qx, M, params.sigma);
    constants[t]->Lr_lut = lut;
    state[t]= new csc_state<T>(M,K);
  }

  int j;
#ifdef _OPENMP
#pragma omp parallel for schedule(static,10),private(j)
#endif
  for (j=0; j < N; j++) {
#ifdef _OPENMP
    int t=omp_get_thread_num();
#else
    int t=0;
#endif
    double* LT = &L[j];
    constants[t]->j = j;
    const csc_constants<T>* constantsT = constants[t];
    csc_state<T>* stateT = state[t];

#if DEBUG
    //    std::cout << "encoding patch " << j << " using thread " << t << std::endl;
#endif
    //
    // for each sample
    //
    // initialize state
    //
    const sc_vector<T> Xj = X.column(j);
    sc_fill(T(0),stateT->a); // a <- all zeroes
    sc_fill(false,stateT->S); // S <- empty
    sc_copy(Xj , stateT->r); // r <- x
    sc_mul(stateT->r,constantsT->D,stateT->Dtr); // Dtr <- D^t*r
    stateT->RSS = sc_dot(stateT->r,stateT->r); 
    stateT->SS = 0;
    stateT->ASA = 0;
    //
    // call core function
    //
    switch(constantsT->step_mode) {
    case CSC_STEP_STEPWISE:
      csc_fss_lg_l_core(constantsT, stateT);
      break;
    case CSC_STEP_MP:
      csc_mp_lg_l_core(constantsT, stateT);
      break;
    case CSC_VANILLA_MP:
      csc_vanilla_mp_lg_l_core(constantsT, stateT);
    }
    //
    // copy results to output
    // 
    *LT = stateT->L; // codelength
    sc_sparse_vector<T> Aj = A.column(j);
    sc_full_to_sparse(stateT->a,Aj);
  }
  //
  // cleanup
  //
#if DEBUG
  //std::cout << "CLEANUP." << std::endl;
#endif
  for (int t=0; t < NTHREADS; t++) {
    delete state[t];    
    delete constants[t]->Lr_lut;
    delete constants[t];
  }
#if DEBUG
  //std::cout << "DONE." << std::endl;
#endif
}
//
//==========================================================================
//==========================================================================
//==========================================================================
//
template<class T>
void csc_fss_lg_l_core(const csc_constants<T>* constants,
		       csc_state<T>* state)
{
  int M = constants->D.m;
  int K = constants->D.n;
  double Lr,Lz,Ls,Lv;
  //
  // sufficient statistics for residual and coefficients
  //
  //
  // initial codelength given by initial solution
  // for this we also keep the empirical residual sum of squares (RSS)
  // 
  double qa = constants->qa;
#ifdef DEBUG
  double qx = constants->qx;
  double L_raw   = M*ceil(log2(2*constants->maxx+1))*qx;
#endif
  // some useful constants
  double sigma2 = constants->sigma * constants->sigma;
  double invM   = 1.0/double(M-1);
  double log2e = log2(exp(1.0));
  double log2M = log2((double)M);
  double log2K = log2((double)K);
  double log2sqrt2PI = log2(sqrt(2.0*M_PI));
  double log2qa = log2(qa);
  std::ofstream a_ostream, state_ostream;

  
  if (constants->dump_path) {
    char fname[256];
    snprintf(fname,255,"A_%06d.ascii",constants->j);
    a_ostream.open(fname);
    snprintf(fname,255,"state_%06d.ascii",constants->j);
    state_ostream.open(fname);
  }

  double theta_r = (invM*state->RSS>sigma2)? sqrt(0.5*(invM*state->RSS-sigma2)) : 0.0;  
  //
  // initial codelength: L(r): L(theta) + L(r|theta)
  //
  Lr = 0.5*log2M + constants->Lr_lut->accu_lookup(state->r.data,state->r.n,theta_r);
  //
  // initial L(a) = describe theta_a = log2K
  //
  Lz = log2K;
  Ls = 0.0;
  Lv = 0.0;
  state->L = Lr + Lz + Ls + Lv;
#ifdef DEBUG
  //std::cout << " Lz=" << log2K << " Ls=" << Ls << " Lv=" << Lv << " L=" << state->L << std::endl;
#endif
  if (constants->dump_path) {
    sc_dump(state->a,a_ostream); // each A is a row 
    state_ostream << state->RSS << ' ' << state->SS << ' ' << Lr << ' ' 
	      << Lz << ' ' << Ls << ' ' << Lv << ' ' << state->L << std::endl;
  }
  //
  // MAIN LOOP: while codelength decreases
  //
  int iter = 0;
#ifdef DEBUG
  //std::cout << std::setprecision(10);
  //std::cout << "ITER=" << iter << " L_raw=" << L_raw << " rss=" << state->RSS << " L=" << state->L;
#endif
  double L_best=1e8;
  double s_best =0.0;
  double a_best =0.0;
  int i_best =-1;
  int SS_best = 0; 
  //double prev_L_best = 1e8;
  double m = 1.0; // qa multiplier for atoms outside active set
  bool stop = false;
  while(!stop) {
#ifdef DEBUGMP
    std::cout << iter;
#endif
    //
    // INNER LOOP: test each candidate
    //
    T* a_data = state->a.data;
    for (int i=0; i < K; i++) {
      double RSS_cand;
      double SS_cand;
      double s_cand = (state->Dtr[i] > 0.0)? qa: -qa;
      //
      // candidate's new i-th coefficient value
      //
      // 
      // step size: allow for bigger step for elements outside active set
      //
      if (!state->S[i])
	s_cand *= m;
      double a_cand = a_data[i] + s_cand;
      //
      // Candidate support size
      //
      if (!state->S[i]) {
	SS_cand = state->SS + 1.0;
      } else if (a_cand == 0.0) { // coef. became 0
	SS_cand = state->SS - 1.0;
      } else {
	SS_cand = state->SS;
      }
      //
      // candidate residual:
      //
      // r^i = r - s^i * Delta * d_k 
      //
      //sc_copy(state->r,state->r_cand);
      sc_scaled_add(1.0,state->r,-s_cand, constants->D.column(i),state->r_cand); 
      //
      // update statistics
      //
      // (r-qa*dk)^T(r-qa*dk) = (r^Tr - 2*qa*r^TD + qa*qa*D^TD)
      // rss(r^i) = rss(r) - 2*s^i*Delta*(D^tr)_i + Delta^2G_{ii}
      //
      RSS_cand = state->RSS - 2.0*s_cand*state->Dtr[i] + s_cand*s_cand*constants->DtD(i,i); 
      // 
      // we are encoding a, conditioned on the fact that |a| >= qa
      // thus, it is enough to encode sgn(a)*max{(|a|-qa),0}
      //
      // PENDING: check whether this can be done more efficiently
      double ac = (fabs(a_data[i]) > qa) ? ((a_data[i]>0.0)? a_data[i]-qa : a_data[i]+qa ) : 0.0;
      double ac_cand = (fabs(a_cand) > qa) ? ((a_cand >0.0)? a_cand-qa : a_cand+qa ) : 0.0;
      // rsa(a^i) = rsa(a) + |ac_i + s_k * Delta| - |ac_i|
      //
      double ASA_cand = state->ASA + fabs(ac_cand) - fabs(ac);
      //
      // compute ML parameters
      //
      theta_r = (invM*RSS_cand > sigma2)? sqrt(0.5*(invM*RSS_cand - sigma2)) : 0.0;  
      //
      // L(r)
      //
      Lr = 0.5*log2M + constants->Lr_lut->accu_lookup(state->r_cand.data, state->r_cand.n, theta_r);
      //
      // L(a): 
      //
      // i) support ~enumerative code: log(K) bits to describe |S|, and log(K choose |S|) for S | |S|
      //    approximated using stirling gives 1/sqrt(2\pi) n^(n+1/2)/( (n-k)^(n-k+1/2) k^(k+1/2) )
      // for S_cand < 3 we do it manually, since the approx is not very good and it is very
      // easy to compute
      if (SS_cand == 1) {
	Lz = log2K + log2K; // log2 (n choose 1) = log2(n)
      } else if (SS_cand == 2) {
	Lz = log2K + log2K + log2(K-1) - log2(2.0); //(n choose 2) = n!/((n-2)!2!) = n*(n-1)/2
      } else if (SS_cand == 3) {
	Lz = log2K + log2K + log2(K-1) + log2(K-2.0) - log2(6.0); //(n choose 3) = n!/((n-3)!3!) = n*(n-1)(n-2)/6
      } else { // Stirling's approximation
	//  with n=K and k=S_cand gives:
	Lz = log2K - log2sqrt2PI + (K + 0.5)*log2(K) - (K-SS_cand+0.5)*log2(K-SS_cand) - (SS_cand+0.5)*log2(SS_cand);
      }
      //
      // ii) sign 
      //
      Ls = SS_cand;
      //
      // iii) magnitude
      // Laplacian, for which an exact formula in terms of rsa exists
      // plus quantization Delta^|S| = -log2(S_cand * log2 qa)
      //
      // PENDING: for large quantization, better encode coefv <- sgn(coefv)*(|coefv|-qa)
      //
      // since the sign has already been taken into account, we need to encode a ONE-sided
      // exponential: L(v) = 1/2 log2(n) -\log P(A|\hat\theta) 
      //                   ~ 1/2 log2(n) -\log p(a|\hat\theta)                        - n\log2\Delta
      //                   = 1/2 log2(n) -\log (1/\hat\theta)^ne^{-\hat\theta||a||_1} - n\log2\Delta
      //                   = 1/2 log2(n) -\log (n/||a||_1)^ne^{-(n/||a||_1)||a||_1}   - n\log2\Delta
      //                   = 1/2 log2(n) + n\log(||a||_1) -n\log n  + n\log e         - n\log2\Delta
      // with n = S_cand 
      //
      // note that because we thresholded a, even if a>0, we may have that the thresholded a
      // used for purposes of coding is all zeroes, resulting in a delta at 0, which is perfectly
      // fine, but gives -Inf when taking log :)
      if (ASA_cand > 0.0)
	Lv = (0.5+SS_cand)*log2(SS_cand) + SS_cand*( log2(ASA_cand) + log2e ) - SS_cand*log2qa;
      else 
	Lv = 0.5*log2(SS_cand); 
      //
      // total
      //
      double L_cand = Lr + Lz + Ls + Lv;
#ifdef DEBUG
      //            std::cout << std::setw(3) << i << ": s=" << std::setw(5) << s_cand << " rss=" << RSS_cand << " rsa=" << ASA_cand << " theta_r=" << theta_r << " Lr=" << Lr << " Lz=" << Lz  << " Ls=" << Ls  
      //      		<< " Lv=" << Lv  << " L=" << L_cand << std::endl;
#endif
      if (L_cand < L_best) { // this is the best i so far
	L_best = L_cand;
	i_best = i;
	s_best = s_cand;
	a_best = a_cand;
	SS_best = SS_cand;
      }	
    } // inner loop
    //
    // best step improves on previous codelength
    //
    if (L_best < state->L) {
      state->L = L_best;
      T prev_a = a_data[i_best];
      a_data[i_best] = a_best; // coefficient
      state->S[i_best] = (a_best != 0.0); // support
      //      for (int j = 0; j < state->S.n; j++) std::cout << ((state->S[j]) ? '1': '0');
      //std::cout << std::endl;
      sc_scaled_add(-s_best, constants->D.column(i_best),state->r);
      sc_vector<T> DtDi = constants->DtD.column(i_best);
      state->RSS += - 2.0*s_best*state->Dtr[i_best] + s_best*s_best*constants->DtD(i_best,i_best); 
      sc_scaled_add(-s_best, DtDi, state->Dtr); 
      state->ASA += (fabs(a_best) - fabs(prev_a));      
      state->SS = SS_best;
      iter++;
#ifdef DEBUG
      //      std::cout << "ITER=" << iter << " L_raw=" << L_raw << " L=" << state->L << " L_best=" << L_best  
      //	      << " i_best=" << i_best << " s_best=" << s_best << " dL= " << (L_best-state->L) << std::endl;
#endif
      if (constants->dump_path) {
	sc_dump(state->a,a_ostream); // each A is a row 
    state_ostream << state->RSS << ' ' << state->SS << ' ' << Lr << ' ' 
	      << Lz << ' ' << Ls << ' ' << Lv << ' ' << state->L << std::endl;
      }
      // reset multiplier to 1
      m = 1.0;
    } else { 
#if 0
      //
      // stopping condition:
      //
      // if by increasing the step size we cannot improve
      // the codelength, relative to the best obtained with the
      // previous smaller step size, then we stop.
      //
      if (m == 1.0) {  // first failure we do not have reference, continue
	prev_L_best = L_best;
	m += 1.0; 
      } else if (L_best < prev_L_best) { 
	//
	// here we increased m and obtained a relative gain
	// keep trying with bigger step: maybe we can reach a point
	// where we improve the current codelength.
	//
	prev_L_best = L_best;
	m += 1.0;
#else
	m += 1.0;
	if (m <= 16.0) {
	  // continue
#endif
      } else {
	//
	// no, we only make things worse by increasing step. STOP
	//
#ifdef DEBUG
	  //	std::cout << "STOP: no further improvement." << std::endl;
#endif
	stop = true;
      }
    }
  } // outer loop
  //
  // clean up/close
  //
  if (constants->dump_path) {
    state_ostream.close();
    a_ostream.close();
  }
}

//
//==========================================================================
//==========================================================================
//==========================================================================
//
template<class T>
void csc_mp_lg_l_core(const csc_constants<T>* constants,
		      csc_state<T>* state)
{
  int M = constants->D.m;
  int K = constants->D.n;
  double sigma2 = constants->sigma * constants->sigma;
  double invM   = 1.0/double(M-1);
  //
  // sufficient statistics for residual and coefficients
  //
  //
  // initial codelength given by initial solution
  // for this we also keep the empirical residual sum of squares (RSS)
  // 
  double qa = constants->qa;
  // components of the codelength
  double Lr,Lz,Ls,Lv;
#ifdef DEBUG
  double qx = constants->qx;
  double L_raw   = M*ceil(log2(2*constants->maxx+1))*qx;
#endif
  // some useful constants
  double log2e = log2(exp(1.0));
  double log2M = log2((double)M);
  double log2K = log2((double)K);
  double log2sqrt2PI = log2(sqrt(2.0*M_PI));
  double log2qa = log2(qa);
  std::ofstream a_ostream, state_ostream;

  if (constants->dump_path) {
    char fname[256];
    snprintf(fname,255,"A_%06d.ascii",constants->j);
    a_ostream.open(fname);
    snprintf(fname,255,"state_%06d.ascii",constants->j);
    state_ostream.open(fname);
  }

  double theta_r = (invM*state->RSS>sigma2)? sqrt(0.5*(invM*state->RSS-sigma2)) : 0.0;  
  //
  // initial codelength: L(r): L(theta) + L(r|theta)
  //
  Lr = 0.5*log2M + constants->Lr_lut->accu_lookup(state->r.data,state->r.n,theta_r);
  //
  // initial L(a) = describe theta_a = log2K
  //
  Lz = log2K;
  Ls = 0.0;
  Lv = 0.0;
  state->L = Lr + Lz;
  
#ifdef DEBUG
  //    std::cout << "Lr= " << Lr << " Lz=" << Lz << " Ls=" << Ls << " Lv=" << Lv << " L=" << state->L << std::endl;
#endif
  if (constants->dump_path) {
    sc_dump(state->a,a_ostream); // each A is a row 
    state_ostream << state->RSS << ' ' << state->SS << ' ' << Lr << ' ' 
	      << Lz << ' ' << Ls << ' ' << Lv << ' ' << state->L << std::endl;
  }
  //
  // MAIN LOOP: while codelength decreases
  //
  //double prev_L_best = 1e8;
  int iter = 0;
  bool stop = false;
  T* a_data = state->a.data;

#ifdef DEBUG
  //  std::cout << std::setprecision(10);
  //  std::cout << "ITER=" << iter << " L_raw=" << L_raw << " rss=" << state->RSS << " L=" << state->L << std::endl;
#endif
  while(!stop) {
#ifdef DEBUGMP
    std::cout << iter;
#endif
  double L_best=1e8;
  double s_best =0.0;
  double a_best =0.0;
  int i_best =-1;
  int SS_best = 0; 
    //
    // INNER LOOP: test each candidate
    //
    for (int i=0; i < K; i++) {
      // does it make sense for it not to be greedy?
      if (state->S[i])
	continue;

      // CSC-MP: step is full correlation, quantized to qa precision
      double s_cand = floor(state->Dtr[i]/qa + 0.5)*qa; 
      //      std::cout << "Dtr[i]=" << state->Dtr[i] << " qa=" << qa << " s_cand=" << s_cand << std::endl;
      if (s_cand == 0.0) {
	// not enough correlation to take a step
	continue;
      }
      //
      // candidate's new i-th coefficient value
      //
      // 
      // step size: allow for bigger step for elements outside active set
      //
      double a_cand = a_data[i] + s_cand;
      //
      // Candidate support size: greedy, always grows!
      //
      //if (!state->S[i]) {
      double SS_cand = state->SS + 1.0;
	//} else if (a_cand == 0.0) { // coef. became 0
	//SS_cand = state->SS - 1.0;
	//} else {
	//SS_cand = state->SS;
	//}
      //
      // candidate residual:
      //
      // r^i = r - s^i * Delta * d_k 
      //
      sc_copy(state->r,state->r_cand);
      sc_scaled_add(-s_cand, constants->D.column(i), state->r_cand); 
      //
      // update statistics
      //
      // (r-qa*dk)^T(r-qa*dk) = (r^Tr - 2*qa*r^TD + qa*qa*D^TD)
      // rss(r^i) = rss(r) - 2*s^i*Delta*(D^tr)_i + Delta^2G_{ii}
      //
      double RSS_cand = state->RSS - 2.0*s_cand*state->Dtr[i] + s_cand*s_cand*constants->DtD(i,i);
      // 
      // we are encoding a, conditioned on the fact that |a| >= qa
      // thus, it is enough to encode sgn(a)*max{(|a|-qa),0}
      //
      // PENDING: check whether this can be done more efficiently
      double ac = (fabs(a_data[i]) > qa) ? ((a_data[i]>0.0)? a_data[i]-qa : a_data[i]+qa ) : 0.0;
      double ac_cand = (fabs(a_cand) > qa) ? ((a_cand >0.0)? a_cand-qa : a_cand+qa ) : 0.0;
      // rsa(a^i) = rsa(a) + |ac_i + s_k * Delta| - |ac_i|
      //
      double ASA_cand = state->ASA + fabs(ac_cand) - fabs(ac);
      //
      // compute ML parameters
      //
      theta_r = (invM*RSS_cand > sigma2)? sqrt(0.5*(invM*RSS_cand - sigma2)) : 0.0;  
      //
      // L(r)
      //
      Lr = 0.5*log2M + constants->Lr_lut->accu_lookup(state->r_cand.data, state->r_cand.n, theta_r);
      //
      // L(a): 
      //
      // i) support ~enumerative code: log(K) bits to describe |S|, and log(K choose |S|) for S | |S|
      //    approximated using stirling gives 1/sqrt(2\pi) n^(n+1/2)/( (n-k)^(n-k+1/2) k^(k+1/2) )
      // for S_cand < 3 we do it manually, since the approx is not very good and it is very
      // easy to compute
#if 1
      if (SS_cand == 1) {
	Lz = log2K + log2K; // log2 (n choose 1) = log2(n)
      } else if (SS_cand == 2) {
	Lz = log2K + log2K + log2(K-1) - log2(2.0); //(n choose 2) = n!/((n-2)!2!) = n*(n-1)/2
      } else if (SS_cand == 3) {
	Lz = log2K + log2K + log2(K-1) + log2(K-2.0) - log2(6.0); //(n choose 3) = n!/((n-3)!3!) = n*(n-1)(n-2)/6
      } else { // Stirling's approximation
	//  with n=K and k=S_cand gives:
	Lz = log2K - log2sqrt2PI + (K + 0.5)*log2(K) 
	  - (K-SS_cand+0.5)*log2(K-SS_cand) - (SS_cand+0.5)*log2(SS_cand);
      }
#else
      // HACK: do not count supp size
      Lz = log2K;
#endif
      //
      // ii) sign 
      //
      Ls = SS_cand;
      //
      // iii) magnitude
      // Laplacian, for which an exact formula in terms of rsa exists
      // plus quantization Delta^|S| = -log2(S_cand * log2 qa)
      //
      // PENDING: for large quantization, better encode coefv <- sgn(coefv)*(|coefv|-qa)
      //
      // since the sign has already been taken into account, we need to encode a ONE-sided
      // exponential: L(v) = 1/2 log2(n) -\log P(A|\hat\theta) 
      //                   ~ 1/2 log2(n) -\log p(a|\hat\theta)                        - n\log2\Delta
      //                   = 1/2 log2(n) -\log (1/\hat\theta)^ne^{-\hat\theta||a||_1} - n\log2\Delta
      //                   = 1/2 log2(n) -\log (n/||a||_1)^ne^{-(n/||a||_1)||a||_1}   - n\log2\Delta
      //                   = 1/2 log2(n) + n\log(||a||_1) -n\log n  + n\log e         - n\log2\Delta
      // with n = S_cand 
      //
      // note that because we thresholded a, even if a>0, we may have that the thresholded a
      // used for purposes of coding is all zeroes, resulting in a delta at 0, which is perfectly
      // fine, but gives -Inf when taking log :)
      if (ASA_cand > 0.0)
	Lv = (0.5+SS_cand)*log2(SS_cand) + SS_cand*( log2(ASA_cand) + log2e ) - SS_cand*log2qa;
      else 
	Lv = 0.5*log2(SS_cand); 
      //
      // total
      //
      double L_cand = Lr + Lz + Ls + Lv;
#ifdef DEBUG
      //      std::cout << std::setw(3) << i << ": s=" << std::setw(5) << s_cand << " rss=" << RSS_cand << " rsa=" << ASA_cand << " theta_r=" << theta_r << " Lr=" << Lr << " Lz=" << Lz  << " Ls=" << Ls  
      //		<< " Lv=" << Lv  << " L=" << L_cand << std::endl;
#endif
      if (L_cand < L_best) { // this is the best i so far
	L_best = L_cand;
	i_best = i;
	s_best = s_cand;
	a_best = a_cand;
	SS_best = SS_cand;
      }	
    } // inner loop
    //
    // best step improves on previous codelength
    //
    iter++;
#ifdef DEBUG
    //    std::cout << "ITER=" << iter << " L_raw=" << L_raw << " L=" << state->L << " L_best=" << L_best  
    //	      << " i_best=" << i_best << " s_best=" << s_best << " dL= " << (L_best-state->L) << std::endl;
#endif
    if (L_best < state->L) {
      double prev_a = a_data[i_best];
      a_data[i_best] = a_best; // coefficient
      state->S[i_best] = (a_best != 0.0); // support
      //      for (int j = 0; j < state->S.n; j++) std::cout << ((state->S[j]) ? '1': '0');
      //std::cout << std::endl;
      sc_scaled_add(-s_best, constants->D.column(i_best), state->r);
      sc_vector<T> DtDi = constants->DtD.column(i_best);
      state->RSS += - 2.0*s_best*state->Dtr[i_best] + s_best*s_best*constants->DtD(i_best,i_best); 
      sc_scaled_add(-s_best, DtDi, state->Dtr); 
      state->ASA += (fabs(a_best) - fabs(prev_a));      
      state->SS = SS_best;
      state->L = L_best;
      if (constants->dump_path) {
	sc_dump(state->a,a_ostream); // each A is a row 
	state_ostream << state->RSS << ' ' << state->SS << ' ' << Lr << ' ' 
		  << Lz << ' ' << Ls << ' ' << Lv << ' ' << state->L << std::endl;
      }
    } else { 
#ifdef DEBUG
      //	std::cout << "STOP: no further improvement." << std::endl;
#endif
	stop = true;
    }
  } // outer loop
  //
  // clean up/close
  //
  if (constants->dump_path) {
    state_ostream.close();
    a_ostream.close();
  }
}


//
//==========================================================================
//==========================================================================
//==========================================================================
// for testing purposes: plian MP
//
template<class T>
void csc_vanilla_mp_lg_l_core(const csc_constants<T>* constants,
		      csc_state<T>* state)
{
  int M = constants->D.m;
  int K = constants->D.n;
  double sigma2 = constants->sigma * constants->sigma;
  double invM   = 1.0/double(M-1);
  //
  // sufficient statistics for residual and coefficients
  //
  //
  // initial codelength given by initial solution
  // for this we also keep the empirical residual sum of squares (RSS)
  // 
  double qa = constants->qa;
  // components of the codelength
  double Lr,Lz,Ls,Lv;
#ifdef DEBUG
  double qx = constants->qx;
  double L_raw   = M*ceil(log2(2*constants->maxx+1))*qx;
#endif
  // some useful constants
  double log2e = log2(exp(1.0));
  double log2M = log2((double)M);
  double log2K = log2((double)K);
  double log2sqrt2PI = log2(sqrt(2.0*M_PI));
  double log2qa = log2(qa);
  std::ofstream a_ostream, state_ostream;
  if (constants->dump_path) {
    char fname[256];
    snprintf(fname,255,"A_%06d.ascii",constants->j);
    a_ostream.open(fname);
    snprintf(fname,255,"state_%06d.ascii",constants->j);
    state_ostream.open(fname);
  }

  double theta_r = (invM*state->RSS>sigma2)? sqrt(0.5*(invM*state->RSS-sigma2)) : 0.0;  
  //
  // initial codelength: L(r): L(theta) + L(r|theta)
  //
  Lr = 0.5*log2M + constants->Lr_lut->accu_lookup(state->r.data,state->r.n,theta_r);
  //
  // initial L(a) = describe theta_a = log2K
  //
  Lz = log2K;
  Ls = 0.0;
  Lv = 0.0;
  state->L = Lr + Lz;
  
#ifdef DEBUG
  //    std::cout << "Lr= " << Lr << " Lz=" << Lz << " Ls=" << Ls << " Lv=" << Lv << " L=" << state->L << std::endl;
#endif
  if (constants->dump_path) {
    sc_dump(state->a,a_ostream); // each A is a row 
    state_ostream << state->RSS << ' ' << state->SS << ' ' << Lr << ' ' 
	      << Lz << ' ' << Ls << ' ' << Lv << ' ' << state->L << std::endl;
  }
  //
  // MAIN LOOP: while codelength decreases
  //
  //double prev_L_best = 1e8;
  int iter = 0;
  bool stop = false;
  T* a_data = state->a.data;

#ifdef DEBUG
  //  std::cout << std::setprecision(10);
  //  std::cout << "ITER=" << iter << " L_raw=" << L_raw << " rss=" << state->RSS << " L=" << state->L << std::endl;
#endif

  while(!stop) {
#ifdef DEBUGMP
    std::cout << iter;
#endif
    //
    // vanilla MP: choose by most correlated, only twist is quantized step
    //
    int i_best = sc_amax(state->Dtr); // O(K)
    T s_best = floor(state->Dtr[i_best]/qa + 0.5)*qa;  // O(1)
    //
    // coefficient value update
    //
    a_data[i_best] += s_best; // O(1)
    //
    // residual update O(M) and correlation update O(K)
    // we wouldn't ever need to update state->r if it wasn't for looking up the codelength!
    //
    sc_scaled_add(-s_best, constants->D.column(i_best), state->r); // O(M)
    state->RSS += constants->DtD(i_best,i_best)*s_best*s_best -2*s_best*state->Dtr[i_best];  // O(1)
    sc_vector<T> DtDi = constants->DtD.column(i_best);
    sc_scaled_add(-s_best, DtDi, state->Dtr);  // O(K)
    //
    // support size update
    //
    state->SS++;
    double SS = state->SS; // shortcut: it is used everywhere
    // 
    // residual codelength
    //
    theta_r = (invM*state->RSS > sigma2)? sqrt(0.5*(invM*state->RSS - sigma2)) : 0.0; // O(1)
    Lr = 0.5*log2M + constants->Lr_lut->accu_lookup(state->r.data, state->r.n, theta_r); // O(M) 

    //
    // coef zeros update and codelength: enumerative code O(1)
    //
    switch (int(SS)) {
    case 1:
	Lz = log2K + log2K; // log2 (n choose 1) = log2(n)
	break;
    case 2:
	Lz = log2K + log2K + log2(K-1) - log2(2.0); //(n choose 2) = n!/((n-2)!2!) = n*(n-1)/2
	break;
    case 3:
      Lz = log2K + log2K + log2(K-1) + log2(K-2.0) - log2(6.0); //(n choose 3) = n!/((n-3)!3!) = n*(n-1)(n-2)/6
      break;
    default:
// Stirling's approximation
	//  with n=K and k=S_cand gives:
	Lz = log2K - log2sqrt2PI + (K + 0.5)*log2(K) 
	  - (K-SS+0.5)*log2(K-SS) - (SS+0.5)*log2(SS);
    }
    //
    // coef sign codelength: O(1)
    //
    Ls = SS;
    //
    // coef magnitude: O(1)
    //
    state->ASA += fabs(s_best);
    //
    // Laplacian, for which an exact formula in terms of rsa exists
    // plus quantization Delta^|S| = -log2(S_cand * log2 qa)
    //
    // PENDING: for large quantization, better encode coefv <- sgn(coefv)*(|coefv|-qa)
    //
    // since the sign has already been taken into account, we need to encode a ONE-sided
    // exponential: L(v) = 1/2 log2(n) -\log P(A|\hat\theta) 
    //                   ~ 1/2 log2(n) -\log p(a|\hat\theta)                        - n\log2\Delta
    //                   = 1/2 log2(n) -\log (1/\hat\theta)^ne^{-\hat\theta||a||_1} - n\log2\Delta
    //                   = 1/2 log2(n) -\log (n/||a||_1)^ne^{-(n/||a||_1)||a||_1}   - n\log2\Delta
    //                   = 1/2 log2(n) + n\log(||a||_1) -n\log n  + n\log e         - n\log2\Delta
    // with n = S_cand 
    //
    // note that because we thresholded a, even if a>0, we may have that the thresholded a
    // used for purposes of coding is all zeroes, resulting in a delta at 0, which is perfectly
    // fine, but gives -Inf when taking log :)
    double effective_ASA = state->ASA - SS*qa; // substract '|qa|' from every coef
    if (effective_ASA > 0.0)
      Lv = (0.5+SS)*log2(SS) + SS*( log2(effective_ASA) + log2e ) - SS*log2qa;
    else 
      Lv = 0.5*log2(SS); 
    //
    // total
    //
    double L_best = Lr + Lz + Ls + Lv;
    iter++;
#ifdef DEBUG
    //    std::cout << "ITER=" << iter << " L_raw=" << L_raw << " L=" << state->L << " L_best=" << L_best  
    //	      << " i_best=" << i_best << " s_best=" << s_best << " dL= " << (L_best-state->L) << std::endl;
#endif
    if (constants->dump_path) {
      sc_dump(state->a,a_ostream); // each A is a row 
      state_ostream << state->RSS << ' ' << SS << ' ' << Lr << ' ' 
		<< Lz << ' ' << Ls << ' ' << Lv << ' ' << state->L << std::endl;
    }
    if (L_best >= state->L) {
#ifdef DEBUG
      //	std::cout << "STOP: no further improvement." << std::endl;
#endif
	stop = true;
    }
  } // outer loop
  //
  // clean up/close
  //
  if (constants->dump_path) {
    state_ostream.close();
    a_ostream.close();
  }
}

#endif
