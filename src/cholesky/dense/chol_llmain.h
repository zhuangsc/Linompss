#ifndef __CHOL_LLMAIN_H__
#define __CHOL_LLMAIN_H__


#ifdef DOUBLE_PRECISION
#define CHOL_HLL			dchol_hll
#define CHOL_LL				dchol_ll
#else
#define CHOL_HLL			schol_hll
#define CHOL_LL				schol_ll
#endif


/* for hypermatrices */
int schol_hll(int mt, int b, int t, float **Ah);  
int dchol_hll(int mt, int b, int t, double **Ah);  


int schol_ll(int mt, int b, int t, float *A, int lda);  
int dchol_ll(int mt, int b, int t, double *A, int lda);  


#endif // __CHOL_LLMAIN_H__
