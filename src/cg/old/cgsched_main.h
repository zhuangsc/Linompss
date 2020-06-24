#ifndef __CGSCHED_MAIN_H__
#define __CGSCHED_MAIN_H__


#include "fptype.h"


#if SINGLE_PRECISION
#define CGSCHED		scgsched
#else
#define CGSCHED		dcgsched
#endif


#define LIBCG_EXPORT __attribute__((__visibility__("default")))


//int       dcg(int bm, int bn, int n, void *A, int s, double *b, double *x, 		     double *xstar, double tol, int steps, double *work, unsigned long works);  
//int       scg(int bm, int bn, int n, void *A, int s, float  *b, float  *x, 			 float  *xstar, float  tol, int steps, float  *work, unsigned long works);  
extern LIBCG_EXPORT int dcgsched(int bm, int bn, int n, void *A, int s, double *b, double *x, int *offs, double *xstar, double tol, int steps, double *work, unsigned long works, int lookahead, int async, double profile, int warmup); 
extern LIBCG_EXPORT int scgsched(int bm, int bn, int n, void *A, int s, float  *b, float  *x, int *offs, float  *xstar, float  tol, int steps, float  *work, unsigned long works, int lookahead, int async, float profile, int warmup);


#endif // __CGSCHED_MAIN_H__

