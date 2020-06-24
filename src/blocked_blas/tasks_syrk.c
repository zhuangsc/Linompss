#include "tasks_syrk.h"

#include "fptype.h"
#include "fpblas.h"
#include "blas.h"


#ifdef SINGLE_PRECISION		

#define __t_syrk		task_ssyrk

#else

#define __t_syrk		task_dsyrk

#endif


void __t_syrk(ompssblas_t uplo, ompssblas_t transa, int n, int k, fp_t alpha, fp_t *A, int lda, fp_t beta, fp_t *C, int ldc, int p) {
	BLAS_syrk(uplo, transa, n, k, alpha, A, lda, beta, C, ldc);
}
