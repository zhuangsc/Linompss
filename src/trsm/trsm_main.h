#ifndef __TRIANGSOLV_MAIN_H__
#define __TRIANGSOLV_MAIN_H__

#include "fptype.h"
#include "fpblas.h"

#ifdef SINGLE_PRECISION
#define TRSM_MAIN	strsm_main
#endif 
#ifdef DOUBLE_PRECISION
#define TRSM_MAIN	dtrsm_main
#endif

int strsm_main(ompssblas_t side, ompssblas_t uplo, ompssblas_t trans, ompssblas_t diag, int m, int n, int b, float alpha, float *A, int lda, float *B, int ldb);

int dtrsm_main(ompssblas_t side, ompssblas_t uplo, ompssblas_t trans, ompssblas_t diag, int m, int n, int b, double alpha, double *A, int lda, double *B, int ldb);

#endif // __TRIANGSOLV_MAIN_H__
