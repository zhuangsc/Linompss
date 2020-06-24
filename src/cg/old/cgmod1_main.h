#ifndef __CGMOD1_MAIN_H__
#define __CGMOD1_MAIN_H__


#if SINGLE_PRECISION
#define CGMOD1		scgmod1
#else
#define CGMOD1		dcgmod1
#endif


int scgmod1(int bm, int bn, int n, void *A, int s, void *b, void *x, int *offs, double tol, int steps, void *work, unsigned long works, int lookahead, int async, double prof, int warmup, int cglog, int release);  
int dcgmod1(int bm, int bn, int n, void *A, int s, void *b, void *x, int *offs, double tol, int steps, void *work, unsigned long works, int lookahead, int async, double prof, int warmup, int cglog, int release);  


#endif // __CGMOD1_MAIN_H__
