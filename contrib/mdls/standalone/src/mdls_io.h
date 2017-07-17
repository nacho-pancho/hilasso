#ifndef MDLS_IO
#define MDLS_IO
#include <sstream>
#include <cstring>
#include <cerrno>
#include "mdls_linalg_types.h"

#define MAXLEN 65536 // fairly large line! accomodates ASCII files of dictionaries with sizes up to 4000 atoms


template<typename T>  
void mdls_read_matlab_matrix(const char* name, const char* filename, 
			     sc_matrix<T>& mat);

template<typename T>
void mdls_read_ascii_matrix(const char* filename, sc_matrix<T>& mat);

template<typename T>
void mdls_write_matlab_matrix(const sc_matrix<T>& mat, const char* name,
			      const char* filename);

template<typename T>
void mdls_write_matlab_matrix(const sc_vector<T>& vec,const char* name,
			      const char* filename);


//
// some matlab-specific structs and tags
//
enum M_FIELD {
  M_LITTLE_ENDIAN = 0
};

enum P_FIELD {
  P_DOUBLE = 0,
  P_FLOAT,
  P_INT,
  P_SHORT,
  P_UNSIGNED,
  P_BYTE
};

enum T_FIELD {
  T_DENSE = 0,
  T_TEXT = 1,
  T_SPARSE = 2
};

struct TypeField {
  char M,O,P,T;
  unsigned int elSize;

  static const unsigned int elemSizes[6];
  TypeField() {}
  TypeField(unsigned int aM,unsigned int aO,
	    unsigned int aP, unsigned int aT){
    M = aM; O = aO; P = aP; T = aT;
  }
  TypeField(long val) {
    M = val % 1000;
    O = (val-M) % 100;
    P = (val-O-M) % 10;
    T = val-P-O-M;
    //const unsigned int TypeField::elemSizes[6] = { 8, 4, 4, 2, 2, 1 };
    switch (P) {
    case 0:
      elSize = 8;
      break;
    case 1: case 2:
      elSize = 4;
      break;
    case 3: case 4:
      elSize = 2;
      break;
    case 5:
      elSize = 1;
      break;
    }
  }
  unsigned int value() { return M*1000 + O*100 + P*10 +T; }
  unsigned int elemSize() { return elSize; }
};

struct Header {
  long type;
  long mrows;
  long ncols;
  long imagf; // imaginary flag
  long namelen;
  Header() {}
  Header(size_t _rows,size_t _cols, long aNamelen) {

    TypeField lType(M_LITTLE_ENDIAN, 0, P_DOUBLE, T_DENSE);
    type = lType.value();
    mrows = _rows;
    ncols = _cols;
    imagf = 0;
    this->namelen = aNamelen;
  }
};


template<typename T>
void mdls_write_matlab_matrix(const sc_matrix<T>& mat,
			      const char* mat_name,
			      const char* filename = NULL)
{
  FILE* lpFile;
  FILE* txtFile;
  Header lHeader;
  
  //TypeField type(M_LITTLE_ENDIAN,0,P_DOUBLE,T_DENSE);
  double *lpColBuff;
  char* lFilename;
  char *txtFilename;
  char *lMatName;
  if (mat_name == NULL) {
    lMatName =  "matrix";
  }
  else {
    lMatName = (char *) mat_name;
  }
  if (filename == NULL){
    lFilename = new char[strlen(lMatName)+5];
    txtFilename = new char[strlen(lMatName)+15];
    strcpy(lFilename,lMatName);
    strcpy(txtFilename,lMatName);
  }
  else {
    lFilename =new char[strlen(filename)+5];
    txtFilename =new char[strlen(filename)+15];
    strcpy(lFilename,filename);
    strcpy(txtFilename,filename);
  }

  if (strstr(lFilename,".mat") == 0)
    strcat(lFilename,".mat");

  lpFile = fopen(lFilename, "wb");
  if (lpFile == 0) {
    std::cerr << lFilename << ":" << strerror(errno) << std::endl;
    return;
  }
  if ( strstr(lFilename,".mat") != 0) {
    int l = strlen(lFilename);
    txtFilename[l-4] = 0;
    strcat(txtFilename,"_mat.txt");
  }
  txtFile = fopen(txtFilename,"w");
  if (txtFile == 0) {
    std::cerr << txtFilename << ":" << strerror(errno) << std::endl;
    return;
  }

  lHeader = Header(mat.m,mat.n,strlen(lMatName)+1);
  lpColBuff = new double[lHeader.mrows];

  fwrite((void*)&lHeader, sizeof(lHeader),1, lpFile);
  fwrite((void*)mat_name, sizeof(char), lHeader.namelen, lpFile);
  for (int j = 0; j < lHeader.ncols ; j++) {
    for (int i = 0; i < lHeader.mrows ; i++) {
      lpColBuff[i] = mat(i,j);
    }
    fwrite((void*)lpColBuff, sizeof(double), lHeader.mrows, lpFile);
  }
  for ( int i=0; i<lHeader.mrows; i++ ) {
    for ( int j=0; j<lHeader.ncols; j++ ) {
      fprintf(txtFile,"%le ",(double)mat(i,j));
    }
    fprintf(txtFile,"\n");
  }
  fclose(lpFile);
  fclose(txtFile);
  delete[] txtFilename;
  delete[] lFilename;
  delete[] lpColBuff;
}

template<typename T>
void mdls_write_matlab_matrix(const sc_vector<T>& v,
			      const char* name,
			      const char* filename)
{
  sc_matrix<T> a(v.data,v.n,1);
  mdls_write_matlab_matrix(a,name,filename);
}


/** 
 * read matrix from ascii text file 
 */
template<typename T>
void mdls_read_ascii_matrix(const char* filename, sc_matrix<T>& mat)
{
  char line[65536];
  double *tmprows[1024]; // max number of rows is this


  std::ifstream in(filename);
  size_t N = 0;
  size_t M = 0;
  size_t n = 0;
  for (int i = 0; i < 1024; i++)
    tmprows[i]=NULL;

  while (!in.eof()) {
    in.getline(line,65536);
    std::istringstream ss(line);
    n = 0;
    if (ss.eof())
      {
        break;
      }
    tmprows[M] = new double[2048];
    while (!ss.eof())
      {
        ss >> tmprows[M][n];
        if (ss.fail())
	  {
	    break;
	  }
        n++;
      }
    if (n == 0) // no numbers here: stop reading
      break;
    if (M == 0)
      {
        N = n;
      } else {
      if (n != N) {
	std::cerr << "Error reading matrix from " << filename << " at line " << M << " column " << n << ": unconsistent no. of cols." << std::endl;
	return;
      }
    }
    M++;
  }
  //
  // fill the matrix
  //
  mat.data = new T[M*N];
  mat.m = M;
  mat.n = N;
  for (size_t i = 0; i < M; i++)  {
      for (size_t j = 0; j < N; j++)
        mat(i,j) = tmprows[i][j];
      delete[] tmprows[i];
    }
  if (tmprows[M] != NULL)
    delete[] tmprows[M];
  //
  // cleanup
  //
}


//
// instantiation for compatible type
//
template<typename T>
void mdls_read_matlab_matrix(const char* name, 
			     const char* filename,
			     sc_matrix<T>& mat)
{
  FILE* lpFile;
  Header lHeader;
  TypeField lType;
  size_t lReadCount;
  size_t lElemCount;
  char *lTmpName;
  double *lpColBuff;
  char* lFilename;

  lFilename = new char[strlen(filename)+5];
  strcpy(lFilename,filename);
  if (strstr(lFilename,".mat") == 0)
    strcat(lFilename,".mat");

  lpFile = fopen(lFilename,"rb");

  if (lpFile == 0) {
    std::cerr << lFilename << ":" << strerror(errno) << std::endl;
    return;
  }

  lReadCount = fread((void*)&lHeader, sizeof(long), 5, lpFile);
#ifdef DEBUG
  std::cout << "Header:\ntype\trows\tcols\timagf\tnamelen\t" << '\n';
  std::cout << "Header:\n" << lHeader.type << '\t' << lHeader.mrows 
	    <<'\t'<<lHeader.ncols<<'\t'<<lHeader.imagf<<'\t'<<lHeader.namelen;
#endif
  lType = TypeField(lHeader.type);
  lTmpName = new char[lHeader.namelen];

  fread((void*)lTmpName, sizeof(char), lHeader.namelen, lpFile);
  lElemCount = lHeader.mrows * lHeader.ncols;
  mat.m = lHeader.mrows;
  mat.n = lHeader.ncols;
  mat.data = new T[mat.m*mat.n];

  lpColBuff = new double[lHeader.mrows];
  for (int j = 0; j < lHeader.ncols; j++) {
    fread((void*)lpColBuff, sizeof(double), lHeader.mrows, lpFile);
    for (int i = 0; i < lHeader.mrows; i++) {
      mat(i,j) = lpColBuff[i];
    }
  }
  fclose(lpFile);
  delete lTmpName;
  delete lFilename;
  delete lpColBuff;
}


#endif
