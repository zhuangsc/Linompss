#ifndef __OMPSS_CG_H__
#define __OMPSS_CG_H__


#ifdef __cplusplus
extern "C" {
#endif

#include "scsparse.h"

#ifdef LIBOMPSS_BUILDING
#define LIBOMPSS_DLL_EXPORTED __attribute__((__visibility__("default")))
#else
#define LIBOMPSS_DLL_EXPORTED
#endif


#if 0
extern LIBOMPSS_DLL_EXPORTED int ompss_scg(int bm, int bn, int n, void *A, int s, float *b, float *x, int *offs, double tol, int steps, void *work, unsigned long works, int lookahead, int async, double prof, int cglog);  

extern LIBOMPSS_DLL_EXPORTED int ompss_dcg(int bm, int bn, int n, void *A, int s, double *b, double *x, int *offs, double tol, int steps, void *work, unsigned long works, int lookahead, int async, double prof, int cglog);  


extern LIBOMPSS_DLL_EXPORTED int ompss_scgmod1(int bm, int bn, int n, void *A, int s, void *b, void *x, int *offs, double tol, int steps, void *work, unsigned long works, int lookahead, int async, double prof, int warmup, int cglog, int release);  

extern LIBOMPSS_DLL_EXPORTED int ompss_dcgmod1(int bm, int bn, int n, void *A, int s, void *b, void *x, int *offs, double tol, int steps, void *work, unsigned long works, int lookahead, int async, double prof, int warmup, int cglog, int release);  
#endif


extern LIBOMPSS_DLL_EXPORTED int ompss_scgprof(int bm, int bn, int n, void *A, int s, void *b, void *x, int *offs, 
		double tol, int steps, void *work, unsigned long works, int lookahead, int async, double profile, int warmup, 
		int cglog, int release, double orth_fac, float *mat_energy, int is_precond, void *zbuff, cs *Acs, css **S, csn **N, int interval, int corrections); 

extern LIBOMPSS_DLL_EXPORTED int ompss_dcgprof(int bm, int bn, int n, void *A, int s, void *b, void *x, int *offs, double tol, 
		int steps, void *work, unsigned long works, int lookahead, int async, double profile, int warmup, int cglog, int release, 
		double orth_fac, float *mat_energy, int is_precond, void *zbuff, cs *Acs, css **S, csn **N, int interval, int corrections); 


#ifdef __cplusplus
}
#endif


#endif // __OMPSS_CG_H__
