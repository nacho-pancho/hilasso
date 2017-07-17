#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "mdls_linalg.h" // linear algebra
#include "img/mdls_pgm.h"   // image reading/writing
#include "img/mdls_deconstruct.h" // image munching
#include "csc/csc.h"        // sparse coding
#include "util/mdls_io.h"  // load the dictionary

char default_dict[] = "ksvd-dict.ascii";
#define NSAMPLES 5

/**
 * This is for profiling purposes
 * 
 * arguments:
 * 1 ......... image file name [lena.pgm]
 * 2 ......... patch width     [8]
 */
int main(int argc, char* argv[]) 
{
  FILE* imgfile;
  const char *dfname, *ifname;
  sc_matrix<double> I, X;
  sc_vector<double> DC, VAR;
  sc_matrix<double> D;

  //
  // parse arguments: 
  //
  // image file
  //
  std::cout << "argc=" << argc << std::endl;
  if (argc < 2)
    ifname = &("lena.pgm")[0];
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
  std::cout << "Loading dictionary " << dfname << std::endl;
  mdls_read_ascii_matrix(dfname,D);
  int M = D.m;
  int K = D.n;
  int w = (int)sqrt(double(M));
#ifdef DEBUG
  std::cout << "M=" << M << " K=" << K << " w=" << w << std::endl;
#endif

  int ov = 0;
  std::cout << "Loading image " << ifname << std::endl;
  //
  // load image
  //
  mdls_pgm_read(imgfile,I.data,I.m,I.n);
  //
  // load dictionary
  //
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
  // choose a few patches and encode them
  //
  std::cout << "Selecting a few patches: ";  
  srand(12345); // very random and unpredictable seed  
  sc_matrix<double> Xs(new double[M*NSAMPLES],M,NSAMPLES);
  for (int i = 0; i < NSAMPLES; i++) {    
    double var = 0.0;
    int randidx = 0;
    while(var < (0.1*maxvar)) {
      randidx = (int) (X.n*double(rand())/double(RAND_MAX));    
      var = VAR[randidx];
    }
    std::cout << randidx << "(var=" << var << "), ";
    sc_vector<double> Xsi = Xs.column(i);
    sc_vector<double> Xi  = X.column(randidx);
    double dc = DC[randidx];
    sc_scaled_add(-dc, Xi, Xsi);
    //
    // remove DC
    //
  }
  std::cout << std::endl;
  csc_coding_params params;
  sc_matrix<double> A(new double[K*NSAMPLES],K,NSAMPLES);
  sc_vector<double> L(new double[NSAMPLES],NSAMPLES);
  std::cout << "Encoding "  << std::endl;  
  params.dump_path = true;
  params.step_mode = CSC_STEP_STEPWISE; 
  csc_fss_lg_l(Xs, D, params, A, L);    
  //
  // free allocated space
  //
  std::cout << "Cleanup " << argv[1] << std::endl;
  delete[] I.data;
  delete[] D.data;
  delete[] X.data;
  delete[] DC.data;
  delete[] VAR.data;
  delete[] A.data;
  delete[] Xs.data;
}
