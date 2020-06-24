#ifndef __CHOL_RLMAIN_H__
#define __CHOL_RLMAIN_H__


#ifdef DOUBLE_PRECISION
#define CHOL_HRL 		dchol_hrl
#define CHOL_RL 		dchol_rl
#else
#define CHOL_HRL 		schol_hrl
#define CHOL_RL 		schol_rl
#endif

/* hypermatrices */
int schol_hrl(int mt, int b, int t, float **Ah);
int dchol_hrl(int mt, int b, int t, double **Ah);


int schol_rl(int mt, int b, int t, float *A, int lda);
int dchol_rl(int mt, int b, int t, double *A, int lda);


#endif // __CHOL_RLMAIN_H__
