#include "trsm_main.h"

#include "bblas_trsm.h"


int TRSM_MAIN(ompssblas_t side, ompssblas_t uplo, ompssblas_t trans, ompssblas_t diag, int m, int n, int b, fp_t alpha, fp_t *A, int lda, fp_t *B, int ldb) 
{

	if ( TEST_LEFT(side) ) {
		printf("ltrsm\n");
		BBLAS_LTRSM(uplo, trans, diag, b, n, m, n, alpha, A, lda, B, ldb);
	} else {
		printf("rtrsm: not yet implemented\n");
	}

	return 0;
}
