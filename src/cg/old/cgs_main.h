#ifndef __CGS_MAIN_H__
#define __CGS_MAIN_H__


#include "hb.h"


#ifdef SINGLE_PRECISION

#define CGS		scgs

#else

#define CGS		dcgs

#endif


int dcgs(int bm, int bn, hbmat_t *A, int s, double *b, double *x, int *offs, double *xstar, double tol, int steps, double *work, unsigned long works, int lookahead, int async, double prof);  
int scgs(int bm, int bn, hbmat_t *A, int s, float  *b, float  *x, int *offs, float  *xstar, float  tol, int steps, float  *work, unsigned long works, int lookahead, int async, float prof);  


#endif // __CGS_MAIN_H__
