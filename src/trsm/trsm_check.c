#include "trsm_check.h"
#include "densutil.h"

#include <math.h>
#include <stdio.h>


extern int m;
extern int n;
extern int lda;
extern int ldb;
extern fp_t alpha;
extern ompssblas_t uplo; /* lower/upper triangular */
extern ompssblas_t side;
extern ompssblas_t trans;
extern ompssblas_t diag;

fp_t trsm_check(int check, fp_t *A, fp_t *X, fp_t *Xchk) 
{
	if ( check ) {
		BLAS_trsm(side, uplo, trans, diag, m, n, alpha, A, lda, Xchk, ldb);
		fp_t err = MAT_NORMDIFF('i', ldb, n, X, ldb, Xchk, ldb);
		printf("check: residule = %.9e\n", err);
		return err;
	}
	return 0;
}
