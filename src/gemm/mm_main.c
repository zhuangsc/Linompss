#include "mm_main.h"

#include "bblas_gemm.h"

int MM_MAIN(int b, int d, int c, int m, int n, int k, ompssblas_t transa, ompssblas_t transb, fp_t alpha, fp_t *A, int lda, fp_t *B, int ldb, fp_t beta, fp_t *C, int ldc) 
{

#ifdef HYPER
	BBLAS_HYPER_GEMM(TRANSP_A, TRANSP_B, b, d, c, m, n, k, FP_ONE, A, m, B, k, FP_NOUGHT, C, m);
#else
	BBLAS_GEMM(1, transa, transb, b, d, c, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
#endif
	
	return 0;
}
