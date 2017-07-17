#include <time.h>
#include "tg_aux.h"

using namespace std;


typedef double cost_array[NGRID][NGRID][DGRID];
///////////////////////////////////////////////////////////////////////
// Arrays where we store the cost values for different metrics,
// indexed by theta, sigma, d
static cost_array mle, l1, l2;
///////////////////////////////////////////////////////////////////////




double
get_utime()
{

	return (double)clock()/CLOCKS_PER_SEC;
}

void fill_arrays(histogram & h, cost_array& mle, cost_array& l1, cost_array& l2,
				 iterbounds t, iterbounds dd, iterbounds s) {
		 double nlg, l1d, l2d;
		 double theta, d, sigma;
		 // fill array of values
		 for ( int i=0; i<t.n; i++ ) {
			 double min_f = 1e100, min_s, min_d;
			 theta = t.s+i*t.i;
			 for ( int j=0; j<s.n; j++ ) {
				 sigma = s.s+j*s.i;
				 for ( int k=0; k<dd.n; k++ ) {
					 d = dd.s+k*dd.i;
					 h.tg_alldist(theta, d, sigma, nlg, l1d, l2d);
					 mle[i][j][k] = nlg;
					 if ( nlg <= min_f ) {
						 // cout << "new min: " << nlg << " s=" << sigma << " d=" << d << endl;
						 min_f = nlg;
						 min_s = sigma;
						 min_d = d;
					 }
					 l1[i][j][k] = l1d;
					 l2[i][j][k] = l2d;
				 }
			 }
#if 0
			 cout << setw(3) << i << " theta=" << setw(8) << theta 
				  << "  min=" << setw(8) << min_f << " s=" << setw(8) << min_s 
				  << " d="<< setw(6) << min_d << " " 
				  << setw(8) <<  get_utime()-time0 << " s" << endl;
#endif
		 }
}

double find_minimum(cost_array A, iterbounds t, iterbounds dd, iterbounds s,				                                 
											  tdsvec& tds,
											  int verbose = 0) {
     double min_nlg=0.0, min_t=0.0, min_d=0.0, min_s=0.0;
	 double sigma, theta, d;
	 int i, j, k;
	 min_nlg = 1e100;
	 for ( i=0, theta=t.s; i<t.n; i++ ) {
		 theta = t.s+i*t.i;
		 for (j=0, sigma=s.s; j<s.n; j++ ) {
			 sigma = s.s + j*s.i;
			 for ( k = 0, d=dd.s; k<dd.n;  k++ ) { 
				 d = dd.s + k*dd.i;
				 double nlg = A[i][j][k];
				 if ( verbose )
					cout << "theta=" << theta << " sig=" << sigma << " d=" << d << " cost=" << nlg << endl;
				if ( nlg <= min_nlg ) {
					min_nlg = nlg;
					min_t = theta;
					min_s = sigma;
					min_d = d;
				}
			 }
		 } 
	 }
	 tds[0] = min_t;
	 tds[1] = min_d;
	 tds[2] = min_s;
	 return min_nlg;
}

void report_minimum(char *label, int classno, double min_cost, double min_t, double min_d, double min_s) {
	printf("%s %3d: Min= %.8le at theta= %6.3lf d= %7.3lf sigma= %7.3lf\n",
		            label,classno,min_cost,min_t,min_d,min_s);
	fflush(stdout);
}

void report_minimum(char *label, int classno, double min_cost, tdsvec &tds) {
	report_minimum(label, classno, min_cost, tds[0], tds[1], tds[2]);
}

void harmonize(double s0, double& s1, double si, int& sn) {

	if ( sn ) s1 = s0 + si*(sn-1);
	else 
		 sn = (int)floor((s1-s0)/si + 1 + si/10);
}

void save_matrix(int classno, cost_type cost,
	                iterbounds t,iterbounds dd, iterbounds s, double min_d) {

		 
		 int k = int(floor((min_d-dd.s)/dd.i+0.5));
		 char *label=NULL;
		 cost_array *ca=NULL;
		 
		 switch(cost) {
		 case MLE:
			 label = "MLE";
			 ca = &mle;
			 break;
		 case L1:
			 label = "L1";
			 ca = &l1;
			 break;
		 case L2:
			 ca = &l2;
			 label = "L2";
			 break;
		 }

		 char fname[256];
		 FILE *ofp;
		 sprintf(fname,"%s/%04d.txt",label,classno);
		 if ( (ofp=fopen(fname,"w"))==NULL ) {
			 cerr << "Error opening file " << fname << " for writing---skipping" << endl;
		 }
		 else {
			 for ( int i=0; i<t.n; i++ ) {
				 for ( int j=0; j<s.n; j++ ) {
					 fprintf(ofp,"%12.8le ",(*ca)[i][j][k]);
				 }
				 fprintf(ofp,"\n");
			 }
			 fclose(ofp);
		 }
 }

double tg_minimize_cost_grid(cost_type cost, histogram h,
			     iterbounds t, iterbounds dd, 
			     iterbounds s, tdsvec &min_tds, int verbose)
{
		 cost_array *ca;
		 switch(cost) {
			 case MLE:
				 ca = &mle;
				 break;
			 case L1:
				 ca = &l1;
				 break;
			 case L2:
				 ca = &l2;
				 break;
		 default:
		   ca = &l2;
		 }
		 fill_arrays(h, *ca, l1, l2, t, dd, s);  // fill 3D arrays of cost values	 
		 // Find MLE (minimum normalized codelength)
		 double min_cost = find_minimum(*ca, t, dd, s, min_tds, verbose);
		 return min_cost;
}
