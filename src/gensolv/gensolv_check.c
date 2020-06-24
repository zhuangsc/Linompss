#include "gensolv_check.h"

#include "fplapack.h"
#include "densutil.h"
#include "matfprint.h"

#include <math.h>
#include <stdio.h>


fp_t gensolv_check(int check, int m, int n, fp_t *A, fp_t *B, fp_t *Aorig, fp_t *Borig) 
{
	if ( check == 0 ) {
		return 0;
	}
	int *ipiv0 = malloc(n*sizeof(int));
	int info;
	LAPACK_getrf(m, n, Aorig, m, ipiv0, &info);
	LAPACK_laswp(n, Borig, m, I_ONE, n, ipiv0, I_ONE);
	BLAS_trsm(OMPSSBLAS_LEFT, OMPSSBLAS_LOWERTRIANG, OMPSSBLAS_NTRANSP, OMPSSBLAS_DIAGUNIT, m, n, FP_ONE, Aorig, m, Borig, m);
	BLAS_trsm(OMPSSBLAS_LEFT, OMPSSBLAS_UPPERTRIANG, OMPSSBLAS_NTRANSP, OMPSSBLAS_NDIAGUNIT, m, n, FP_ONE, Aorig, m, Borig, m);
	fp_t err = MAT_NORMDIFF('i', m, n, B, m, Borig, m);
	printf("check: residual norm = %.9e\n", err);

	return err;
}
