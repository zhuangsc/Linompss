#ifndef __CGPROF_MAIN_H__
#define __CGPROF_MAIN_H__


#include "fptype.h"
#include "scsparse.h"


#ifdef SINGLE_PRECISION
#define CGPROF		scgprof
#else
#define CGPROF		dcgprof
#endif


#define LIBCG_EXPORT __attribute__((__visibility__("default")))


extern LIBCG_EXPORT int dcgprof(int bm, int bn, int n, void *A, int s, void *b, void *x, int *offs, double tol, 
		int steps, void *work, unsigned long works, int lookahead, int async, double profile, int warmup, int cglog, 
		int release, double orth_fac, float *mat_energy, int is_precond, void *zbuff, cs *Acs, css **S, csn **N, int interval, int corrections); 

extern LIBCG_EXPORT int scgprof(int bm, int bn, int n, void *A, int s, void *b, void *x, int *offs, double tol, 
		int steps, void *work, unsigned long works, int lookahead, int async, double profile, int warmup, int cglog, 
		int release, float orth_fac, float *mat_energy, int is_precond, void *zbuff, cs *Acs, css **S, csn **N, int interval, int corrections);


#endif // __CGPROF_MAIN_H__

