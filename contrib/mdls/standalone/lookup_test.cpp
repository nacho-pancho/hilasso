#include <iostream>
#include <fstream>
#include "src/lookup.h"
#include <cstdlib>

int main(int argc, char* argv[]) {
  int M = 64;
  int maxx = 255;
  double sigma = 2.0; // very little noise
  //double beta = 256.0*0.07;
  double kappa = 2.8;

  if (argc > 1)
    sigma = atof(argv[1]);

  parametric_lut* lut;
#ifdef DEBUG
  std::cout << "initializing Laplacian LUT."  << std::endl;
#endif
  double qx = 1.0; 
  std::ofstream fout;

  lut = create_laplacian_parametric_lut(maxx,qx,M);
  fout.open("lut-lap.ascii");
  lut->dump(fout);
  fout.close();
  delete lut;
  lut = NULL;

#ifdef DEBUG
  std::cout << "initializing MOE LUT."  << std::endl;
#endif
  lut = create_moe_parametric_lut(maxx,qx,M,kappa);
  fout.open("lut-moe.ascii");
  lut->dump(fout);
  fout.close();
  delete lut;
  lut = NULL;

#ifdef DEBUG
  std::cout << "initializing LG LUT."  << std::endl;
#endif
  lut = create_lg_parametric_lut(maxx,qx,M,sigma);
  fout.open("lut-lg.ascii");
  lut->dump(fout);
  fout.close();
  delete lut;
  lut = NULL;

#ifdef DEBUG
  std::cout << "initializing MOEG LUT."  << std::endl;
#endif
  lut = create_moeg_parametric_lut(maxx,qx,M,kappa,sigma);
  fout.open("lut-moeg.ascii");  
  lut->dump(fout);
  fout.close();
  delete lut;
  lut = NULL;

#ifdef DEBUG
  std::cout << "Accessing LUT."  << std::endl;
#endif
  
#ifdef DEBUG
  std::cout << "destroying LUT."  << std::endl;
#endif
  delete lut;
#ifdef DEBUG
  std::cout << "done."  << std::endl;
#endif
}

