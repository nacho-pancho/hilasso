#ifndef MDLS_PGM
#define MDLS_PGM

#include <cstdio>
#include <iostream>

#define MAXROWLEN 8192
#define MAXROWBUF 16384
#define LINELEN MAXROWLEN

typedef unsigned char buffer_t;

//----------------------------------------------------------------------------
// READ
//----------------------------------------------------------------------------

/**
 * I=pgmread(fname)
 *
 * inputs (prhs):
 * 0 fname: image name
 *
 * outputs (plhs):
 * 0 I : image
 */
int read_pnm_header (FILE * infile, int *colsp, int *rowsp, int *maxsp,
		     int *typep);


template <typename T>
void mdls_pgm_read(FILE* imgfile, T*& I, int& M, int& N) 
{
  unsigned char rowbuf[MAXROWBUF];
  int lsb,msb;
  int maxs, type;
  int L = M*N;
  int nr;
  read_pnm_header(imgfile, &N,&M,&maxs,&type);
  std::cout << "M=" << M << " N=" << N << std::endl;

  if ( (I == NULL) || (L==0) ) {
    if (maxs < 256) {
      I = new T[M*N];
    } else {
      I = new T[M*N*2];
    }
  } else {
    if (maxs < 256) {
      if (L < (M*N)) {
	std::cerr << "I is not large enough, reallocating." << std::endl;
	delete[] I;
	I = new T[M*N];
      }
    } else {
      if (L < (2*M*N)) {
	std::cerr << "I is not large enough, reallocating." << std::endl;
	delete[] I;
	I = new T[2*M*N];
      }
    }
  }  
  if (maxs < 256) { /* 16 bit values */
    for (int i = 0; i<M; i++) {
      nr = fread(&rowbuf[0],sizeof(unsigned char),N,imgfile);
      for (int j = 0; j < N; j++)
	I[M*j+i] = T(rowbuf[j]);
    }
  } 
  else { /* 16 bit values */
    for (int i = 0; i<M; i++)
    { 
      nr = fread(&rowbuf[0],sizeof(char),2*N,imgfile);
      for (int j = 0; j < N; j++) {
        /* MSB then LSB */
	msb = rowbuf[2*j];
	lsb = rowbuf[2*j+1];
        I[M*j+i] = T(((msb<<8)&0xff00) | (lsb&0xff)); // col-major access
      }
    }
  }
}

//----------------------------------------------------------------------------
// WRITE
//----------------------------------------------------------------------------
template<typename T>
void mdls_pgm_write(T*& I, int M, int N, FILE* imgfile) 
{
  unsigned char rowbuf[MAXROWBUF];
  double scale;

  int MN   = M*N;
  double minv =  1.0e15;
  double maxv = -1.0e15;
  //
  // scan range of values
  //
  for (int i=0; i < MN; i++) {    
    if (I[i] > maxv)
      maxv=I[i];
    else if (I[i] < minv)
      minv=I[i];
  }
  //
  // compute normalization factor to fit in 0-255 range
  //
  if ((maxv < 256.0) && (minv >=0)) {
    scale=1.0;
    minv=0.0;
  } else {
    scale=255.0/(maxv-minv);
  }
  //
  // write header  in chunks of MAXROWBUF bytes
  //  
  fprintf(imgfile,"P5 %d %d 255\n",M,N);
  size_t k=0;
  for (int i=0; i < M; i++) {
    for (int j=0; j < N; j++) {
      rowbuf[k++] = (unsigned char)((double(I[i+j*M])-minv)*scale +0.5);
      if (k >= MAXROWBUF) {
        fwrite(rowbuf,sizeof(unsigned char),MAXROWBUF,imgfile);
        k=0;
      }
    }
  }
  //
  // write trailing bytes
  //
  if (k > 0) {
    fwrite(rowbuf,sizeof(unsigned char),k,imgfile);
  }
}
#endif
