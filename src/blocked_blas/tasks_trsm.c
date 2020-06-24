#include "tasks_trsm.h"

#include "fptype.h"
#include "fpblas.h"

#include "matfprint.h"


#ifdef SINGLE_PRECISION
#define __t_trsm			task_strsm
#endif

#ifdef DOUBLE_PRECISION
#define __t_trsm			task_dtrsm
#endif


// C = alpha * A * B - beta * C
// A bm x bk 
// B bk x bn 
// C bm x bn 
void __t_trsm(ompssblas_t side, ompssblas_t uplo, ompssblas_t trans, ompssblas_t diag, int bm, int bn, fp_t alpha, fp_t *A, int lda, fp_t *B, int ldb, int p) 
{
	BLAS_trsm(side, uplo, trans, diag, bm, bn, alpha, A, lda, B, ldb);
}

