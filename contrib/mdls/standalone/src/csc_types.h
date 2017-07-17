#ifndef MDLS_CSC_TYPES
#define MDLS_CSC_TYPES

#include "mdls_linalg.h"
#include "lookup.h"

template<class T>
class csc_state { // current state for CSC-FSS algorithm
 public:

  sc_vector<T> a;    // current coefficients
  sc_vector<T> r;    // current residual
  sc_vector<T> Dtr;  // current correlation
  sc_vector<bool> S;    // current support  
  sc_vector<T> r_cand;   // candidate 
  sc_vector<T> Dtr_cand;   // candidate 
  double RSS; // current sum of squares of residual
  double ASA;  // current sum of absolute of coefficients
  double L ; // current codelength
  double SS; // current support size

  csc_state(int M, int K) {
    a.allocate(K);
    r.allocate(M);
    Dtr.allocate(K);
    r_cand.allocate(M);    
    Dtr_cand.allocate(K);
    S.allocate(K);
  }

  ~csc_state() { 
    a.free();
    r.free();
    Dtr.free();
    r_cand.free();
    Dtr_cand.free();
    S.free();
  }
};

struct csc_cd_state { // current state for CSC-CD algorithm
  // pending
};

/**
 * variants of step selection in CSC FSS algorithm
 * the CSC_VANILLA_MP uses only correlation, not codelength,
 * to select a new factor. Faster but gives significantly 
 * larger codelengths. Only for benchmarking.
 */
typedef enum {CSC_STEP_STEPWISE = 0, CSC_STEP_MP = 1, CSC_VANILLA_MP = 2} step_modes;

/**
 * parameters of the CSS family of algorithms
 * must be passed through this structure. This avoids
 * errors due to wrong argument placement in function calls.
 * The parameters are as follows:
 *
 * maxx ....... number of possible values of  samples in X (integer, usually 255)
 * qx ........ quantization step for X. Set to <=0 for default.
 * qa ........       ,,      ,, ,,  A.  ,,
 * sigma ..... std. dev of Gaussian noise. If <=0, the std. dev. of the corresponding 
 *             quantization error, qx/sqrt(12) will be used.
 * step_mode .. stepsize selection mode. See step_modes.
 * dump_path .... if true, files will be generated that store the state of the algorithm
 *            at each step. The files will be called sampleN_<var>.ascii where
 *            <var> will be A for the coefficients, R for the residual, L for codelength and S for
 *            the support size.
 *
 */
struct csc_coding_params {

  double qx; // signal quantization step
  double qa; // coefficients quantization step  
  double maxx;   // x takes values in [0,maxx]
  double sigma; // variance of Gaussian noise in model
  bool dump_path; // if true, large files will be generated with the solution path and intermediate results obtained during coding of each sample. Use with care.
  int step_mode; // 0: non-greedy, stepwise, 1: greedy matching pursuit
  int nzmax; // forced maximum number of nonzero elements in reconstruction
  double Loffset; // additional fixed codelength to be added to the final result

  csc_coding_params() {
    maxx = 255.0;
    qx = 1.0;
    qa = 256.0/16.0;
    sigma = 1.0/sqrt(12.0);
    step_mode = CSC_STEP_STEPWISE;
    dump_path = false;
    Loffset = 0.0;
    nzmax = 0;
  }

};

#endif
