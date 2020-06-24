#ifndef __LU_MAIN_H__
#define __LU_MAIN_H__

#include "fptype.h"

#ifdef SINGLE_PRECISION
#define LU_MAIN		slu_main
#endif

#ifdef DOUBLE_PRECISION
#define LU_MAIN		dlu_main
#endif

void slu_main(int m,int n,int b, float *A, int lda, int *IPIV);
void dlu_main(int m,int n,int b, double *A, int lda, int *IPIV);

#endif // __LU_MAIN_H__
