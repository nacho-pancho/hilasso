#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "src/mdls_linalg.h" // linear algebra
#include "src/mdls_pgm.h"   // image reading/writing
#include "src/mdls_deconstruct.h" // image munching
#include "src/csc.h"        // sparse coding
#include "src/mdls_io.h"  // load the dictionary
#include "src/mdls_dict_learn.h" // learn dictionary

char default_dict[] = "ksvd-dict.ascii";

/**
 * 
 * arguments:
 * 1 ......... image file name [barbara.pgm]
 * 2 ......... dictionary file name [ksvd-dict.ascii]
 * 3 ......... number of iterations [10]
 */
int main(int argc, char* argv[]) 
{
  FILE* imgfile;
  const char *dfname, *ifname;
  sc_matrix<double> I, X;
  sc_vector<double> DC, VAR;
  sc_matrix<double> D;
  std::ofstream out;
  //
  // parse arguments: 
  //
  // image file
  //
  std::cout << "argc=" << argc << std::endl;
  if (argc < 2)
    ifname = &("barbara.pgm")[0];
  else
    ifname = argv[1];
  imgfile = fopen(ifname,"r");
  if ( imgfile == NULL) {
    std::cerr << "Unable to read " << ifname << std::endl;
    exit(1);
  }
  //
  // dictionary file
  //
  if (argc < 3)
    dfname = default_dict;
  else
    dfname = argv[2];
  //
  // number of iterations
  //
  int J = 20;
  if (argc >= 4)
    J = atoi(argv[3]);
  std::cout << "Will perform " << J << " iterations of DL " << std::endl;
  int ov = 0;
  std::cout << "Loading image " << ifname << std::endl;
  //
  // load image
  //
  mdls_pgm_read(imgfile,I.data,I.m,I.n);
  //
  // load dictionary
  //
  std::cout << "Loading dictionary " << dfname << std::endl;
  mdls_read_ascii_matrix(dfname,D);
  int M = D.m;
  int K = D.n;
  int w = (int)sqrt(double(M));
  std::cout << "M=" << M << " K=" << K << " w=" << w << std::endl;
  //
  // decompose image
  //
  std::cout << "Deconstructing into patches "  << std::endl;
  bool computeDC = true;
  bool computeVAR = true;
  mdls_deconstruct(I,w,ov,computeDC,computeVAR,X,DC,VAR);
  double maxvar = 0.0;
  for (int i = 0; i < VAR.n; i++)
    if (VAR[i] > maxvar)
	    maxvar = VAR[i];
  std::cout << "Maximum variance=" << maxvar << std::endl;

  //
  // remove DC and discard constant patches
  // (variance below 1/100 of maximum variance)
  //
  int N = X.n;
  int Ns = 0;
  for (int i = 0; i < N; i++) {    
    if (VAR[i] >= .01*maxvar) Ns++;
  }
  std::cout << "Kept " << Ns << " patches for learning, discarded " << N-Ns << std::endl;
  sc_matrix<double> Xs(new double[M*Ns],M,Ns);
  for (int i = 0, is=0; i < N; i++) {
    if (VAR[i] >= .01*maxvar) {
      sc_vector<double> x  = X.column(i);
      sc_vector<double> xs = Xs.column(is);
      for (int j=0; j < x.n; j++) {
        xs[j] = x[j] - DC[i]; // without DC
      }
      is++;
    }
  }
  out.open("Xs.ascii");
  sc_dump(Xs,out);
  out.close();
  //
  // adapt dictionary to image
  //
  int Ntest = Ns/5;
  int Ntrain  = Ns - Ntest;
  // these are shallow copies of parts of Xs
  sc_matrix<double> Xtrain(&Xs.data[0],M,Ntrain);
  sc_matrix<double> Xtest(&Xs.data[M*Ntrain],M,Ntest); 
  csc_learning_data<double> data(Xtrain,Xtest);
  
  srand(12345); // very random and unpredictable seed  
  csc_learning_params mparams;  
  mparams.max_iter = J;
  mparams.batch_size = 500;
  mparams.debug = DEBUG_MILD;
  mparams.xval_step = 5;
  mparams.dump = true;
  mparams.inertia = 0.8;

  csc_coding_params cparams;
  cparams.nzmax = K;
  cparams.step_mode = CSC_STEP_MP;
  csc_learning_state<double>* pstate = new csc_learning_state<double>(M,K,J);
  // 
  // learn using L2
  //
  std::cout << "Learning using L2 error model" << std::endl;
  if (true) {
    sc_copy(D,pstate->D);
    mdls_dict_learn_csc_l2(data,cparams,mparams,*pstate);
    out.open("csc_l2_dict.ascii");
    sc_dump(pstate->D,out);
    out.close();
  }
  delete pstate;
  //
  // learn using Huber
  //
  std::cout << "Learning using Huber error model" << std::endl;
  pstate = new csc_learning_state<double>(M,K,J);
  sc_copy(D,pstate->D);
  mdls_dict_learn_csc_huber(data,cparams,mparams,*pstate);
  out.open("csc_huber_dict.ascii");
  sc_dump(pstate->D,out);
  out.close();
  //
  // free allocated space
  //
  std::cout << "Cleanup " << argv[1] << std::endl;
  delete[] I.data;
  delete[] D.data;
  delete[] X.data;
  delete[] DC.data;
  delete[] VAR.data;
  delete[] Xs.data;

  delete pstate;
}
