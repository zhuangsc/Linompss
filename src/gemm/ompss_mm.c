#include "ompss_mm.h"

#include "mm_main.h"
//#include "bblas_gemm.h"

int ompss_dgemm(ompssblas_t transa, ompssblas_t transb, int b, int d, int c, int m, int n, int k, double alpha, const double *A, int lda, const double *B, int ldb, double beta, double *C, int ldc) 
{
	dmm_main(b, d, c, m, n, k, transa, transb, alpha, A, lda, B, ldb, beta, C, ldc);
//	bblas_dgemm(1, transa, transb, b, d, c, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

int ompss_sgemm(ompssblas_t transa, ompssblas_t transb, int b, int d, int c, int m, int n, int k, float alpha, const float *A, int lda, const float *B, int ldb, float beta, float *C, int ldc) 
{
	smm_main(b, d, c, m, n, k, transa, transb, alpha, A, lda, B, ldb, beta, C, ldc);
//	bblas_sgemm(1, transa, transb, b, d, c, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}


/* d and c are not supported yet */
int ompss_scsrmb(int b, int d, int c, int n, float alpha, hbmat_t *Acsr, const float *B, int ldb, float beta, float *C, int ldc) 
{
	if ( Acsr->b != b || b == 0 ) {
		return 1;
	}

	scsrmmb(b, d, c, n, alpha, Acsr, B, ldb, beta, C, ldc); 

	return 0;
}

/* d and c are not supported yet */
int ompss_dcsrmb(int b, int d, int c, int n, double alpha, hbmat_t *Acsr, const double *B, int ldb, double beta, double *C, int ldc) 
{
	if ( Acsr->b != b || b == 0 ) {
		return 1;
	}

	dcsrmmb(b, d, c, n, alpha, Acsr, B, ldb, beta, C, ldc); 

	return 0;
}


