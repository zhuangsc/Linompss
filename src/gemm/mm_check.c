#include "matmul_check.h"

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>


int gemm_check(int check, int b, int d, int c, int m, int n, int k, ompssblas_t transa, ompssblas_t transb, fp_t alpha, 
		void *A, int lda, void *B, int ldb, fp_t beta, void *C, int ldc) {
	if ( check ) {

		fp_t *Cx = malloc( ldc * max(m, n) * sizeof(fp_t) );

		BLAS_gemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, Cx, ldc);

		fp_t norm = MAT_NORMDIFF('i', m, n, C, m, Cx, m);

		printf("check: diff norm %e\n", norm);

		free(Cx);
	}

	return 0;
}
