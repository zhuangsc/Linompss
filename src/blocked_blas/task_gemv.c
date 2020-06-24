#include "task_gemv.h"

#include "fpblas.h"
#include "fptype.h"
#include "fpmatr.h"


#if DOUBLE_PRECISION
#define __t_gemv 	task_dgemv
#else
#define __t_gemv 	task_sgemv
#endif


void __t_gemv(int bm, int bn, int m, int n, fp_t alpha, fp_t *A, fp_t *x, fp_t beta, fp_t *y) {
	BLAS_gemv(OMPSSBLAS_NTRANSP, bm, bn, alpha, A, m, x, i_one, beta, y, i_one);
}
