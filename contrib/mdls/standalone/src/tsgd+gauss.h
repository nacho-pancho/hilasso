#ifndef _TSGD_GAUSS_
#define _TSGD_GAUSS_

#include <iostream>
#include <iomanip>
#include <math.h>
#include <gsl/gsl_vector.h>
#include "tg_probs.h"


using namespace std;

#ifndef ABSIZE
	#define ABSIZE 256
#endif

#ifndef MULTIMIN_N
#define MULTIMIN_N 3
#endif

#define MIN_THETA 1e-4
#define MAX_THETA (1.-MIN_THETA)


typedef enum { MLE, L1, L2} cost_type;

 

// bounds for a parameter
class iterbounds {
public:
	double s, e, i; // start, end, increment
	int    n;       // number of points

	iterbounds() { s=e=i=0.0; n=0; }
	iterbounds(double ss, double ee, double ii, int nn=0) {
		s = ss;
		e = ee;
		i = ii;
		n = nn;
		if ( n ) e = s + i*(n-1);
		else 
		 n = (int)floor((e-s)/i + 1 + i/10);
		if ( s>e || i==0 ) {
			cerr << "Bad bounds---setting e=s, n=1" << endl;
			e = s;
			n = 1;
			i = 1;
		}
	}

	void print() {
		cout << "min=" << s << " max=" << e << " incr=" << i << " num=" << n << endl;
	}
	void print(char *label) {
		cout << label << ": ";
		print();
	}
	void enforce(double& x) {
		if ( x<s ) x = s;          // clip to lower bound
		else if ( x > e ) x = e;   // clip to upper bound
	}
};
// vector definition for minimization calculations
class tdsvec {
public:
	double tds[MULTIMIN_N];

	inline tdsvec() { tds[0] = tds[1] = tds[2] = 0.0; }

	inline tdsvec(double a, double b, double c) {
		tds[0] = a; tds[1] = b; tds[2] = c;
	}

	inline void set_theta(double v) { tds[0]=v; }
	inline void set_d(double v) { tds[1]=v; }
	inline void set_sigma(double v) { tds[2]=v; }

	inline double get_theta() const { return tds[0]; }
	inline double get_d() const { return tds[1]; }
	inline double get_sigma() const { return tds[2]; }

	inline double norm() { 
		double f = 0;
		for ( int i=0; i< MULTIMIN_N; i++ ) f += fabs(tds[i]);
		return f;
	}

	double& operator[](int i) { return tds[i];	}

	const double& operator[](int i) const { return tds[i]; }

	void mulscalar(double a) { 
		for ( int i=0; i< MULTIMIN_N; i++ ) tds[i] *= a;
	}
	void add(tdsvec y) {
		for ( int i=0; i< MULTIMIN_N; i++ ) tds[i] += y[i];
	}
	void print() {
		for ( int i=0; i<MULTIMIN_N; i++ ) 
			cout << setw(6) << tds[i] << " ";
	}

	friend std::ostream& operator<<(std::ostream& out, const tdsvec& tds) {
	  out << "theta=" << setw(5) << tds[0] << " d=" << setw(5) << tds[1] << " s=" << setw(5) << tds[2];
	  return out;
	}

	void normalize(double min_sigma) {
		if ( tds[0] < 0 ) tds[0] = 0;
		else if ( tds[0] > 1 ) tds[0] = 1;
		if ( tds[2] < min_sigma ) tds[2] = min_sigma;
	}
	int steepest() {
		double f0 = fabs(tds[0]),
			   f1 = fabs(tds[1]),
			   f2 = fabs(tds[2]);
		if ( f0 >= f1 && f0 >= f2 )
			return 1-2*(tds[0]<0);
		else 
			if ( f1 >= f2 ) 
				return 2 - 4*(tds[1]<0);
			else
				return 3 - 6*(tds[2]<0);
	}
	void from_gsl_vec(gsl_vector *x) {
		for ( int i=0; i< MULTIMIN_N; i++ )
			tds[i] = gsl_vector_get(x, i);
	}
	void to_gsl_vec(gsl_vector *x) const {
		for ( int i=0; i< MULTIMIN_N; i++ )
			gsl_vector_set(x, i, tds[i]);
	}

};



class histogram {
	int hist[2*ABSIZE-1];
	int N;
public:
	void init() {
		for ( int i=0; i<2*ABSIZE-1; i++ ) hist[i]=0;
		N = 0;
	}
	histogram() { init(); }
	histogram(int *pt) {
		int *p = pt;
		init();
		for ( int i=-(ABSIZE-1); i<ABSIZE; i++,p++) {
			(*this)[i] = *p;
			N += *p;
		}
	}
	histogram(histogram &h) {
		N = 0;
		for ( int i=-(ABSIZE-1); i<ABSIZE; i++ ) {
			(*this)[i] = h[i];
			N += h[i];
		}
	}

	int& operator[](int i) {
		return hist[ABSIZE-1+i];
	}

	const int& operator[](int i) const {
		return hist[ABSIZE-1+i];
	}

	int getn() {
		return N;
	}

	double avg() {
		int s = 0;
		int n = N;
		for ( int i=-(ABSIZE-1); i<ABSIZE; i++ ) {
			int x = (*this)[i];
			s += i*x;
		}
		return (double)s/n;
	}

	double var() {
		int s = 0;
		int n = N;
		for ( int i=-(ABSIZE-1); i<ABSIZE; i++ ) {
			s += i*i*(*this)[i];
		}
		double E = avg();
		return (double)s/n-(E/n)*(E/n);
	}

	double avgabs() {
		int s = 0;
		int n = N;
		for ( int i=1; i<ABSIZE; i++ ) {
			s += i*((*this)[i]+(*this)[-i]);
		}
		return (double)s/n;
	}

	void increment(int x, int delta=1) {
		N += delta;
		if ( x <= -ABSIZE || x >= ABSIZE ) {
			cerr << "histogram.increment: index " << x << " out of range---ignored" << endl; 
			return;
		}
		(*this)[x] += delta;
	}



	void print() {
		int n = N;
		for ( int i = -(ABSIZE-1); i < ABSIZE; i++ ) {
			if ( (*this)[i] ) {
				double p = double((*this)[i])/n;
				cout << setw(5) << i << " " <<setw(8) << (*this)[i] << setw(10) << " "<< p << endl;
			}
		}
		cout << "Total samples: " << n << endl;
		if ( n != N ) {
			cout << "*********** INCONSISTENT: N = " << N << " *** ABORTING ***"<< endl;
			exit(-10);
		}
		
	}

	friend inline std::ostream& operator<<(std::ostream& out, const histogram& h) {
	  out << '[';
	  for (size_t i = 0; i < 511; i++)
	    out << setw(5) << h.hist[i] << ", ";
	  out << ']';
	  return out;
	}

	double tg_prob_log(double theta, double d, double sigma)
	{
		double lg = 0.0;

		for ( int i = 1-ABSIZE; i < ABSIZE; i++ ) {
			if ( (*this)[i] ) {
				double p = TSGD_Gauss_P(i, theta, d, sigma);
				lg += (*this)[i]*log(p);
			}
		}
		return lg/log(2.0);
	}

	double tg_norm_codelength(double theta, double d, double sigma) {
		return -tg_prob_log(theta, d, sigma)/N;
	}

	double tg_L2dist(double theta, double d, double sigma) {
		double dist = 0.0;
		for ( int i = 1-ABSIZE; i < ABSIZE; i++ ) {
			if ( (*this)[i] ) {
				double p = TSGD_Gauss_P(i,  theta, d, sigma);
				double dif = p - (double)(*this)[i]/N;
				dist += dif*dif;
			}
		}
		return dist;
	}

	double tg_L1dist(double theta, double d, double sigma = 0) {
		double dist = 0.0;
		for ( int i = 1-ABSIZE; i < ABSIZE; i++ ) {
			if ( (*this)[i] ) {
				double p = TSGD_Gauss_P(i, theta, d, sigma);
				double dif = fabs(p - (double)(*this)[i]/N);
				dist += dif;
			}
		}
		return dist;
	}

	void tg_alldist(double theta, double d, double sigma, 
		double& ncodelength, double&L1, double & L2)
	{

		double lg = 0.0;
		L1 = L2 = 0.0;

		for ( int i = 1-ABSIZE; i < ABSIZE; i++ ) {
			if ( (*this)[i] ) {
				double p = TSGD_Gauss_P(i, theta, d, sigma);
				lg += (*this)[i]*log(p);
				double dif = fabs(p - (double)(*this)[i]/N);
				L1 += dif;
				L2 += dif*dif;
			}
		}
		ncodelength = -lg/(log(2.0)*N);
	}


	double cost_fxn(cost_type cost, double theta, double d, double sigma) {
		switch(cost) {
			case MLE:
				return tg_norm_codelength(theta, d, sigma);
				break;
			case L2:
				return tg_L1dist(theta, d, sigma);
				break;
			case L1:
				return tg_L2dist(theta, d, sigma);
				break;
			default:
				cerr << "Unknown cost specification " << cost << endl;
				exit(10);
		}
		return 0.0;
	}

	double cost_fxn(cost_type cost, tdsvec x) {
		return cost_fxn(cost, x[0], x[1], x[2] );
	}

#define EPS 0.001  // epsilon for gradient computation

	int tg_gradient(cost_type cost, double theta, double d, double sigma, tdsvec& g)
	{
		double eps = EPS;
		double t0, t1;
		// gradient at x is computed by taking differences between the values at x+eps and x-eps.
		// theta
		t0 = theta-eps; if (t0 <= MIN_THETA ) t0 = MIN_THETA;
		t1 = theta+eps; if ( t1 >= MAX_THETA ) t1 = MAX_THETA;
		g[0] = (cost_fxn(cost,t1, d, sigma) - cost_fxn(cost,t0, d, sigma))/(t1-t0);
		
					
		// d (not bounded)
		g[1] = (cost_fxn(cost,theta, d+eps, sigma) - cost_fxn(cost,theta, d-eps, sigma))/(2*eps);

		//sigma
		t0 = sigma-eps; 
		t1 = sigma+eps;
		g[2] = (cost_fxn(cost,theta, d, t1) - cost_fxn(cost,theta, d, t0))/(t1-t0);

		return 0;
	}

	int tg_gradient(cost_type cost, tdsvec& tds, tdsvec& g) {
		return tg_gradient(cost, tds[0], tds[1], tds[2], g);
	}

#if 0   // not in use, obsolete
	double tg_minimize_cost_descent(cost_type cost, tdsvec &tdsinit, tdsvec& tdsmin,  double min_sigma,
		                                             double step=1e0, double tol=1e-5) {

			 tdsvec tds = tdsinit, 
				    tds1,
				    g;
			 

			 tds.normalize(min_sigma);
			 tg_gradient(cost, tds, min_sigma, g);
			 
			 int iter = 0;
			 double nrm0 = g.norm(), 
				    nrm1, 
					eps = step;
			 while( nrm0 > tol ) {

				 iter++;

				 g.mulscalar(-eps);
				 tds.add(g);
				 tds.normalize(min_sigma);
				 tg_gradient(cost, tds, min_sigma, g);
				 nrm1 = g.norm();
				 double ng2 = nrm1-nrm0;
				 cout << setw(5) << iter; 
				 cout << "  g: "; g.print(); cout << " N: " << setw(6) << nrm1 << " " << setw(6) << ng2/nrm1 << endl;
				 nrm0 = nrm1;
				 
			 }
			 return cost_fxn(cost, tds);
	}


	double tg_minimize_cost_greedy(cost_type cost, tdsvec &tdsinit, tdsvec& tdsmin,  double min_sigma,
		                                             tdsvec step) {
			 tdsvec x = tdsinit;
			 double fmin;
			 int done = 0, iter=0;
			 while ( !done ) {
				 tdsvec g;
				 iter++;
				 tg_gradient(cost, x, min_sigma, g);
				 int k = g.steepest();
				 fmin = cost_fxn(cost, x);
				 cout << setw(5) << iter; 
				 cout << "  x: "; x.print(); cout << " k: " << setw(2) <<  k <<
					 " f: " << setw(6) << fmin << endl;
				 if ( k<0 ) x[-k-1] += step[-k-1];
				 else x[k-1] -= step[k-1];
			 }
			 return cost_fxn(cost, x);
	}


#endif


	int read(FILE *fp) {
		init();
		for ( int j=-(ABSIZE-1); j<ABSIZE; j++ ) {
			int ct;
			double dct;
			int rc = fscanf(fp, "%le",&dct);
			if ( rc != 1 || feof(fp) ) { 
				if ( j > -(ABSIZE-1) ) // OK to fail on first entry
				cerr << "histogram.read:**** Failed on entry " << j << endl;
				return -1;
			}
			// counts file can contain floating point or integers
			ct = (int)floor(dct+0.5); // convert to int
			increment(j, ct);
		}
		return 0;
	}
	
};

#define MAXROWS (1<<16)

class histogram_array {
	histogram *mx[MAXROWS];
	int ncols;
	int nrows;
public:
	void init() { ncols=nrows=0; }
	histogram_array(){ init();}
	int read(char *fname, int nc=2*ABSIZE-1) 
	{
		FILE *fp;
		//		int totct = 0;
		ncols = nc;
		init();
		if ( (fp=fopen(fname,"r"))==NULL ) {
			perror(fname);
			return -1;
		}
		while ( 1 ) {
			histogram h;
			if ( h.read(fp) < 0 ) {
				break;
			}
			mx[nrows] = new histogram[1];
			*mx[nrows] = h;
			nrows++;
			if ( nrows >= MAXROWS ) {
				cerr << "histogram_array.read: max number of rows reached " << nrows << endl;
				break;
			}
		}
		return nrows;
	}
	histogram& operator[](int i) { 
		if ( i<0 || i>=nrows ) {
			cerr << "histogram_array: index " << i << " out of range 0.." << nrows-1 << endl;
			throw(-1);
		}
		return *mx[i];
	}
	int get_nrows() { return nrows; }
};

#endif

