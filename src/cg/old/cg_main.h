#ifndef __CG_MAIN_H__
#define __CG_MAIN_H__


#ifdef SINGLE_PRECISION

#ifdef USE_SPARSE
#define CG		bcn_scgs
#else
#define CG		bcn_scg
#endif

#else

#ifdef USE_SPARSE
#define CG		bcn_dcgs
#else
#define CG		bcn_dcg
#endif

#endif




int bcn_dcg(int bm, int bn, int n, void *A, int s, double *b, double *x, int *offs, double tol, int steps, double *work, unsigned long works, int lookahead, int async, double prof, int cglog);  
int bcn_dcgs(int bm, int bn, int n, void *A, int s, double *b, double *x, int *offs, double tol, int steps, double *work, unsigned long works, int lookahead, int async, double prof, int cglog);  
int bcn_scgs(int bm, int bn, int n, void *A, int s, float  *b, float  *x, int *offs, float  tol, int steps, float  *work, unsigned long works, int lookahead, int async, float prof, int cglog);  
int bcn_scg(int bm, int bn, int n, void *A, int s, float  *b, float  *x, int *offs, float  tol, int steps, float  *work, unsigned long works, int lookahead, int async, float prof, int cglog);  


#endif // __CG_MAIN_H__

