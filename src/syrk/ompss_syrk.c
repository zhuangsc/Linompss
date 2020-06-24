#include "ompss_syrk.h"
#include "syrk_main.h"

void ompss_ssyrk(ompssblas_t uplo, ompssblas_t trans, int b, int n, int k, float alpha, float *A, int lda, float beta, float *C, int ldc)
{
	ssyrk_main(uplo, trans, b, n, k, alpha, A, lda, beta, C, ldc);
}

void ompss_dsyrk(ompssblas_t uplo, ompssblas_t trans, int b, int n, int k, double alpha, double *A, int lda, double beta, double *C, int ldc)
{
	dsyrk_main(uplo, trans, b, n, k, alpha, A, lda, beta, C, ldc);
}

