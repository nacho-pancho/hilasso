#ifndef _TG_AUX_
#define _TG_AUX_

#include "tsgd+gauss.h"
#define NGRID  200
#define DGRID  22

void save_matrix(int classno, cost_type cost, 
				 iterbounds t,iterbounds dd, iterbounds s, double min_d);
void harmonize(double s0, double& s1, double si, int& sn);
void harmonize(double s0, double& s1, double si, int& sn);
void report_minimum(char *label, int classno, double min_cost, tdsvec &tds);
void report_minimum(char *label, int classno, double min_cost, double min_t, double min_d, double min_s);
double tg_minimize_cost_grid(cost_type cost, histogram h, 
								   iterbounds t, iterbounds dd, iterbounds s, 
								   tdsvec &min_tds, int verbose=0);
#endif
