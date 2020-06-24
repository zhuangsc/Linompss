#ifndef __MATMUL_CHECK_H__
#define __MATMUL_CHECK_H__

#include "fptype.h"
#include "fpblas.h"
#include "fpmatr.h"
#include "densutil.h"

int gemm_check(int check, int b, int d, int c, int m, int n, int k, ompssblas_t transa, ompssblas_t transb, fp_t alpha, void *A, int lda, void *B, int ldb, fp_t beta, void *C, int ldc); 

int matmul_check(int check, int b, int d, int c, int m, int n, int k, int transa, int transb, void *A, void *B, int ldb, void *C, int ldc);
#endif // __MATMUL_CHECK_H__
