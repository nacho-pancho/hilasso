#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include "mdls_linalg.h"
#include <gsl/gsl_interp.h>

#ifdef DEBUG
#define my_assert(a) assert((a))
#else
#define my_assert(a) 
#endif

/**
 * Parametric density function.
 */
struct sc_function {
  double (*fun)(double x, double* params);
  double* params;
};
/**
 * Lookup table for uniformly quantized variables (for fast indexing) on a
 * bounded symmetric interval and precomputed derivatives. Derivatives are
 * computed using GSL cubic interpolation functions.
 * This is intended for symmetric probability density functins  f(x) on the real
 * line and related functions such as -log f(x) and [-log f(x)]'
 * The function is supplied by the user.
 */
class lut {
 
public:
 lut(void): values(NULL),derivatives(NULL),tsize(0) { }

  /**
   * initializes the table with user supplied function 
   */ 
  lut(double _xamax, double dx, sc_function& pdf) {
    idx = 1.0/dx;
    tsize = 2*int(ceil(xmax*idx))+1;
    offset = int(ceil(xmax*idx));
    xmax = double(offset)*dx;
    xmin = -xmax;
    values = new double[tsize];
    derivatives = new double[tsize];
    double x = -xmax;
    for (int i = 0; i < tsize; i++, x += dx) {
      values[i] = (pdf.fun)(x,pdf.params);
    }
    _compute_derivative();
  }

  /**
   * initializes table with custom data. User must provide 2*ceil(xamax/dx)+1 values
   */
  lut(double xamax, double dx, double* _values) { 
    idx = 1.0/dx;
    tsize = 2*int(ceil(xmax*idx))+1;
    offset = int(ceil(xmax*idx));
    xmax = double(offset)*dx;
    xmin = -xmax;
    values = new double[tsize];
    derivatives = new double[tsize];
    for (int i = 0; i < tsize; i++) {
      values[i] = _values[i];
    }
    _compute_derivative();
  }
  
  ~lut() { 
    if (values!=NULL) {
      delete[] values; 
      delete[] derivatives;
    }    
  }
  
  inline double operator[](int i) const {
    return values[i];
  }
  
  inline double& operator[](int i)  {
    return values[i];
  }
  
  inline double val(int i) const {
    return values[i+offset];
  }
  
  inline double val(double x) const {
    return values[int(x*idx+0.5)+offset];
  }
  
  inline double der(int i) const {
    return derivatives[i+offset];
  }
  
  inline double der(double x) {
    return derivatives[int(x*idx+0.5)+offset];
  }
  
  //  private:
  void _compute_derivative() {
    double *x;
    gsl_interp* p_interp_info = NULL;
    gsl_interp_accel *accel = gsl_interp_accel_alloc();
    
    p_interp_info = gsl_interp_alloc(gsl_interp_cspline ,tsize);
    x = new double[tsize];
    double dx = 1.0/idx;
    double xi = xmin;
    for (int i = 0; i < tsize; i++, xi+= dx) {
      x[i] = xi;
    }
    gsl_interp_init(p_interp_info,x,values,tsize);
    xi = xmin;
    for (int i = 0; i < tsize; i++, xi+= dx) {
      derivatives[i] = gsl_interp_eval_deriv(p_interp_info,x,values,xi,accel);
    }
    gsl_interp_free(p_interp_info);
    gsl_interp_accel_free(accel);
    delete[] x;
  }

  double* values; /// values of the function
  double* derivatives; /// derivatives of the function
  int tsize; /// table size
  int offset; /// offset of x=0
  double idx,xmax,xmin; /// inverse of dx, max and min values represented
};

/**
 * Lookup table class for encoding two-part universal probability models
 * where the parameter is quantized in 
 * step if 1/sqrt(M), where M is the length of the sequence (x_1,x_2,...,x_M) to encode.
 * And the x_i take values in an alphabet [-maxv:qx:maxv ]
 *
 * TODO: use lut class for each line
 */
class parametric_lut {
 public:

  /**
   * Initialize tables for parametric distributions of samples x which take values in a
   * discrete support [-maxx:qx:maxx]. Thus each table has size 2*maxx/qx+1 
   * The parameters for which the tables are instantiated take values from [0:qt:maxt],
   * for a total of [maxt/qt+1] tables.
   */
  parametric_lut(double _maxx, double _qx, double _maxt, double _qt) {
    maxx = _maxx;
    minx = maxx;
    maxt = _maxt;
    iqx = 1.0/_qx;
    offset = (int)ceil(maxx*iqx);
    tsize = 2*offset+1;
    iqt = 1.0/_qt;
    ntables = ceil(maxt*iqt) + 1;
    //ntables = ceil(sqrt(maxx*double(M)));
    data = new double*[ntables];
    for (int i=0; i < ntables; i++) {
      data[i] = new double[tsize];
    }    
  }

  /**
   * Initialize from externally created table. Useful for tables created in Matlab or similar.
   * The first row of the external table has the support of x.
   * The first column has the values of theta corresponding to the table on each row.
   */
  parametric_lut(const sc_matrix<double>& ext_lut) {
    sc_vector<double> v;
    tsize = ext_lut.n-1; // first column does not count
    ntables = ext_lut.m-1; // first row does not count
    //
    // retrieve LUT parameters from 'headers'
    //
    v = ext_lut.column(0);
    maxt = v[ntables]; // last row
    iqt  = 1.0/(maxt-v[ntables-1]);
    v    = ext_lut.column(tsize); // last column
    maxx = v[0];  // -maxx:
    v    = ext_lut.column(tsize-1); // next to last column
    iqx  = 1.0/(maxx - v[0]);
    offset = (int)ceil(maxx*iqx);


    //
    // copy data
    //
    data = new double*[ntables];
    for(int i=0; i < ntables; i++) {
      data[i] = new double[tsize];
      for(int j=0; j < tsize; j++) {
	data[i][j] = ext_lut(i+1,j+1); // skip first column and row
      }
    }
    //
    // done
    //
  }

  ~parametric_lut() {
    for (int i = 0; i < ntables; i++) {
      delete[] data[i];
    }
    delete[] data;
  }

  inline double lookup(double x,double theta) {
    my_assert( (x >= minx) && (x <= maxx));
    return data[ get_table_index(theta) ][ offset + (int)(iqx*x) ];
  }

  void dump(std::ofstream& fout) {
    // first row is x values
    // first col is theta values
    double qx = 1.0/iqx;
    fout << "0 ";
    for (int j = 0; j < tsize; j++) {
      fout << double(j)*qx-double(offset) << ' ';
    }
    fout << '\n';
    for (int i = 0; i < ntables; i++) {
      fout << double(i)/iqt << ' ';
      for (int j = 0; j < tsize; j++) {
	fout << data[i][j] << ' ';
      }
      fout << '\n';
    }
  }

  /** 
   * batch evaluation: replaces each element in x by f(x;theta)
   */
  inline void batch_lookup(double* x, int n, double theta) {
    double* tmptable = data [get_table_index(theta)];
    for (int i = 0; i < n; i++)
      x[i] = tmptable[ offset + (int)(iqx*x[i]) ];
  }

  /** 
   * cumulative evaluation: result = sum_i f(x_i;theta)
   */
  inline double accu_lookup(const double* x, int n, double theta) {
    double* tmptable = data [get_table_index(theta)];
    double y = 0.0;
    for (int i = 0; i < n; i++) {
      int j = offset + (int)(iqx*x[i]);
      if (j < 0)
	j = 0;
      else if (j >= tsize)
	j = tsize-1;
      y += tmptable[ j ];
    }
    return y;
  }

  inline double get_qt() { return 1.0/iqt; }
  inline double get_qx() { return 1.0/iqx;  }
  inline double** get_data() const { return data; } // warning: raw access
  inline int get_ntables() const { return ntables;  }
  inline int get_tsize() const { return tsize;  }
  inline int get_offset() const { return offset;  }
  inline double get_maxx() const  { return maxx; }
  inline double get_minx() const { return minx; }
  inline double get_maxt() const { return maxt; }

 private:

  double **data;
  double iqt; // table i indexes table for parameter value theta=qt*i
  double iqx; // table index j indexes value x=qx*j-offset
  double maxt;//maximum parameter value (minimum assumed 0)
  double minx, maxx;
  int offset; // offset of x=0 in table
  int ntables;
  int tsize;

  inline int get_table_index(double theta) {
    int ti = (int)(iqt*theta);
    return (ti < ntables) ? ti : (ntables-1);
  }

};


/**
 * create lookup tables for the discretized Laplacian model
 *
 * maxx .......... maximum positive value of samples.                  
 * qx ............ Quantization step of samples.
 * M ............. number of samples to encode. This gives the quantization
 *                 step for the parameter theta
 */
parametric_lut* create_laplacian_parametric_lut(int maxx, double qx, int M);

/**
 * create lookup tables for the discretized MOE model
 *
 * maxx .......... maximum positive value of samples.                  
 * qx ............ Quantization step of samples.
 * M ............. number of samples to encode. This gives the quantization
 *                 step for the parameter theta
 * kappa ......... MOE shape parameter.
 */
parametric_lut* create_moe_parametric_lut(int maxx, double qx, int M, double kappa);

/**
 * create lookup tables for the discretized Laplacian+Gaussian 
 * model (LG). 
 *
 * maxx .......... maximum positive value of samples.                  
 * qx ............ Quantization step of samples.
 * M ............. number of samples to encode. This gives the quantization
 *                 step for the parameter theta
 * sigma ......... Gaussian std. dev.
 */
parametric_lut* create_lg_parametric_lut(int maxx, double qx, int M,double sigma2);

/**
 * create lookup tables for the discretized Laplacian+Gaussian 
 * model (LG). 
 *
 * maxx .......... maximum positive value of samples.                  
 * qx ............ Quantization step of samples.
 * M ............. number of samples to encode. This gives the quantization
 *                 step for the parameter theta
 * kappa ......... MOE shape parameter
 * sigma ......... Gaussian std. dev.
 */
parametric_lut* create_moeg_parametric_lut(int maxx, double qx, int M, double kappa, double sigma);

#undef my_assert
