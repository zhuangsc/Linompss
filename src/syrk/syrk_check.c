#include "syrk_check.h"

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

fp_t syrk_check(int check, ompssblas_t uplo, ompssblas_t trans, int b, int n, int k, fp_t alpha, fp_t *A, int lda, fp_t beta, fp_t *C, int ldc, fp_t *Cc) 
{
	if ( check ) {

		BLAS_syrk(uplo, trans, n, k, alpha, A, lda, beta, Cc, ldc);

		fp_t rerr = DMAT_RELERR(uplo, n, n, C, Cc);

		printf("check: rel err %e \n", rerr);

		return rerr;

	}

	return 0;
}
