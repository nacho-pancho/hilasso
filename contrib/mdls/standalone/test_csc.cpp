#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "src/mdls_linalg.h" // linear algebra
#include "src/mdls_pgm.h"   // image reading/writing
#include "src/mdls_deconstruct.h" // image munching
#include "src/mdls_reconstruct.h" // image munching
#include "src/csc.h"        // sparse coding
#include "src/mdls_io.h"  // load the dictionary
#include <gsl/gsl_cblas.h>

char default_dict[] = "ksvd-dict.ascii";
#define NSAMPLES 20

/**
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
  time_t t0,t1;

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
  fclose(imgfile);
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
  int N = X.n;
  double maxvar = 0.0;
  for (int i = 0; i < VAR.n; i++)
    if (VAR[i] > maxvar)
	    maxvar = VAR[i];
  std::cout << "Maximum variance=" << maxvar << std::endl;
  //
  // remove dc
  //
  for (int j = 0; j < N; j++) {
    for (int i = 0; i < M; i++) {
      X.data[M*j+i] -= DC[j];
    }
  }
  //
  // for analysis of the algorithm: 
  // choose a few patches and encode them, dump info to files
  //
  std::cout << "Encoding a few patches for analysis ";  
  srand(12345); // very random and unpredictable seed  
  sc_matrix<double> Xs;
  Xs.allocate(M,NSAMPLES);

  for (int i = 0; i < NSAMPLES; i++) {    
    double var = 0.0;
    int randidx = 0;
    while(var < (0.1*maxvar)) {
      randidx = (int) (N*double(rand())/double(RAND_MAX));    
      var = VAR[randidx];
    }
    std::cout << randidx << "(var=" << var << "), ";
    sc_vector<double> Xsi = Xs.column(i);
    sc_vector<double> Xi  = X.column(randidx);
    sc_copy(Xi,Xsi);
  }
  std::cout << std::endl;
  //
  // initialize coding parameters
  //
  csc_coding_params params;
  params.nzmax = M/4;
  std::cout << "nzmax=" << params.nzmax << std::endl;
  params.dump_path = true;
  // 
  // allocate structures for result
  //
  sc_sparse_matrix<double> A;
  sc_vector<double> L;
  A.allocate(K,NSAMPLES,params.nzmax);
  L.allocate(NSAMPLES);
  std::cout << "Encoding "  << std::endl;  

  csc_fss_lg_l(Xs, D, params, A, L);    
  std::cout << "Overall codelength: " << sc_sum(L) << std::endl;
  std::ofstream dump_ostream;

  dump_ostream.open("Xs.ascii");
  sc_dump(Xs,dump_ostream);
  dump_ostream.close();
  
  dump_ostream.open("As.ascii");
  sc_dump(A,dump_ostream);
  dump_ostream.close();

  // recovered: D*A
  sc_scaled_mac(1.0,D,false,A,false,0.0,Xs);
  dump_ostream.open("Xr.ascii");
  sc_dump(Xs,dump_ostream);
  dump_ostream.close();

#if 0
  //
  // big multiplication: D^TX, done in GSL and my functions to see the difference
  //
  sc_matrix<double> DtX;
  DtX.allocate(K,N);
  time(&t0);
  for (int k=0; k < 10; k++) {
    sc_mul(D,true,X,false,DtX);
  }
  time(&t1);
  std::cout << "sc_mul time " << t1-t0 << std::endl;
  //
  // CBLAS way
  //
  time(&t0);
  for (int k=0; k < 10; k++) {   
    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,
		K,M,M,1.0,D.data,M,X.data,M,0.0,
		DtX.data,M);
  }
  time(&t1);
  std::cout << "cblas time " << t1-t0 << std::endl;
  delete[] DtX.data;
#endif
  //
  // for profiling: disable dump, encode all patches
  //
#if 1
  A.free();
  L.free();
  A.allocate(K,N,params.nzmax);
  L.allocate(N);
  sc_matrix<double> Xr,Ir;
  Xr.allocate(X.m,X.n);
  Ir.allocate(I.m,I.n);

  params.dump_path = false;
  for (int sm = 0; sm <= 2; sm++) {
    char aux[100];
    std::cout << "Encoding all " << N << " patches for analysis using mode " << sm << std::endl;  
    time(&t0);
    params.step_mode = sm;
    csc_fss_lg_l(X, D, params, A, L);    
    
    time(&t1);
    // compressed length: what csc returns plus describing the DC, which
    // takes log2(256*w^2) per wxw block, and there are N of those
    double CL = sc_sum(L) + log2(256.0*w*w)*N;
    double UL = double(N)*w*w*8.0; // raw uncompressed length
    int nnz = 0;
    for (int j=0; j < A.n; j++) {
      sc_sparse_vector<double> Aj = A.column(j);
      nnz += Aj.nnz();      
    }
    std::cout << "Average nnz=" << ((double)nnz)/((double)A.n) << std::endl;
    std::cout << "Overall codelength: " << sc_sum(L) << std::endl;
    std::cout << "Compression ratio: " << CL/UL << std::endl;
    double dt = t1-t0;
    std::cout << "done in " << dt << "s. Average: " << double(N)/double(dt) << " patches/s." << std::endl;
    snprintf(aux,99,"A_mode%d.ascii",sm);
    dump_ostream.open(aux);
    sc_dump(A,dump_ostream);
    dump_ostream.close();
    // 
    // compose approximate sparse image (Excluding residual)
    //
    if (sm==1) {
      sc_scaled_mac(1.0,D,false,A,false,0.0,Xr);
      mdls_reconstruct(Xr.data,w,Xr.n,
		       Ir.m,Ir.n,ov,
		       DC.data,Ir.data,(double*)NULL);
      
      imgfile = fopen("approx.pgm","w");
      mdls_pgm_write(Ir.data,Ir.m,Ir.n,imgfile);
      fclose(imgfile);
    }
    snprintf(aux,99,"L_mode%d.ascii",sm);
    dump_ostream.open(aux);
    sc_dump(L,dump_ostream);
    dump_ostream.close();
  }
  Xr.free();
#endif
  //
  // free allocated space
  //
  std::cout << "Cleanup " << argv[1] << std::endl;
  I.free();
  D.free();
  X.free();
  DC.free();
  VAR.free();
  A.free();
  L.free();
  Xs.free();
}
