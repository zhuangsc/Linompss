#ifndef __MM_H__
#define __MM_H__

#include "fptype.h"
#include "fpmatr.h"

#ifdef DOUBLE_PRECISION

#define MM_MAIN		dmm_main

#endif

#ifdef SINGLE_PRECISION

#define MM_MAIN		smm_main

#endif


int dmm_main(int b, int d, int c, int m, int n, int k, ompssblas_t transa, ompssblas_t transb, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc);

int smm_main(int b, int d, int c, int m, int n, int k, ompssblas_t transa, ompssblas_t transb, float  alpha, float  *A, int lda, float  *B, int ldb, float  beta, float  *C, int ldc);


#endif // __MM_H__
