#include "chol_check.h"

#include <math.h>
#include <stdlib.h>

#include "fptype.h"
#include "fplapack.h"
#include "fpblas.h"
#include "densutil.h"


int chol_check(int check, int m, int mr, int ts, fp_t **Lh, fp_t *L, fp_t *A) 
{
	if( check > 0 ) {
		printf( "check: potrf...");
		fflush(0);

		fp_t *L1 = malloc(m*m*sizeof(fp_t));
		fp_t *L2 = malloc(m*m*sizeof(fp_t));
		LAPACK_lacpy(LAPACK_TRIANG_LOWER, m, m, L, m, L1, m);
		LAPACK_lacpy(LAPACK_TRIANG_LOWER, m, m, L, m, L2, m);
		BLAS_trmm(OMPSSBLAS_RIGHT, OMPSSBLAS_LOWERTRIANG, OMPSSBLAS_TRANSP, OMPSSBLAS_NDIAGUNIT, m, m, FP_ONE, L1, m, L2, m);
#if 0
		int info;
		LAPACK_potrf(OMPSSLAPACK_LOWERTRIANG, m, A, m, info);
		printf( "...done (info: %d)\n", info);
		if ( info < 0 ) {
			fprintf(stderr, "check: potrf error\n");
			return 1;
		}
#endif

		fp_t norm_diff = DMAT_RELERR(OMPSSBLAS_LOWERTRIANG, m, m, L2, A);
		printf( "norm of difference: %22.16g\n", norm_diff);

		int ret = (norm_diff > 0.0001) ? 1 : 0 ;

		return ret;
	}

	return 0;
}
