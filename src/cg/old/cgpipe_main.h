#ifndef __CGPIPE_MAIN_H__
#define __CGPIPE_MAIN_H__


#if SINGLE_PRECISION

#if USE_SPARSE
#define CGPIPE			scgpipes
#else
#define CGPIPE			scgpipe
#endif

#else

#if USE_SPARSE
#define CGPIPE			dcgpipes
#else
#define CGPIPE			dcgpipe
#endif

#endif


/* sparse */
int scgpipes(int bm, int bn, int n, void *A, int s, float *b, float *xbuff, int *offs, float tol, int steps, float *work, unsigned long works, int lookahead, int async, float profile, int warmup, int cglog); 
int dcgpipes(int bm, int bn, int n, void *A, int s, double *b, double *xbuff, int *offs, double tol, int steps, double *work, unsigned long works, int lookahead, int async, double profile, int warmup, int cglog); 
/* dense */
int scgpipe(int bm, int bn, int n, void *A, int s, float *b, float *xbuff, int *offs, float tol, int steps, float *work, unsigned long works, int lookahead, int async, float profile, int warmup, int cglog); 
int dcgpipe(int bm, int bn, int n, void *A, int s, double *b, double *xbuff, int *offs, double tol, int steps, double *work, unsigned long works, int lookahead, int async, double profile, int warmup, int cglog); 


#endif // __CGPIPE_MAIN_H__
