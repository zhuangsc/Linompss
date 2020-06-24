#include "mm.h"

#include "bblas_gemm.h"
#include "fptype.h"


int MM(int p, int b, int d, int c, int m, int n, int k, ompssblas_t transa, ompssblas_t transb, fp_t *A, fp_t *B, fp_t *C) {
#ifdef HYPER
	BBLAS_HYPER_GEMM(b, d, c, m, n, k, FP_ONE, A, m, B, k, FP_NOUGHT, C, m);
#else
	BBLAS_GEMM(p, OMPSSBLAS_NOTRANSP, OMPSSBLAS_NOTRANSP, b, d, c, m, n, k, FP_ONE, A, m, B, k, FP_NOUGHT, C, m);
#endif
	
	return 0;
}
