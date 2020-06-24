#include "ompss_trsm.h"
#include "trsm_main.h"

void ompss_strsm(ompssblas_t side, ompssblas_t uplo, ompssblas_t trans, ompssblas_t diag, int m, int n, int b, float alpha, float *A, int lda, float *B, int ldb)
{
	strsm_main(side, uplo, trans, diag, m, n, b, alpha, A, lda, B, ldb);
}

void ompss_dtrsm(ompssblas_t side, ompssblas_t uplo, ompssblas_t trans, ompssblas_t diag, int m, int n, int b, double alpha, double *A, int lda, double *B, int ldb)
{
	dtrsm_main(side, uplo, trans, diag, m, n, b, alpha, A, lda, B, ldb);
}
