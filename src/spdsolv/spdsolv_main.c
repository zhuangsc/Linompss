#include "spdsolv_main.h"

#include "fptype.h"
#include "fpblas.h"
#include "fpmatr.h"
#include "bblas_trsm.h"
#include "ompss_dense_chol.h"


int SPDSOLV_MAIN(int m, int n, int b, fp_t *A, fp_t *B) 
{

	OMPSS_CHOL_RL(m, b, b, A, m);

	BBLAS_LTRSM(OMPSSBLAS_LOWERTRIANG, OMPSSBLAS_NTRANSP, OMPSSBLAS_NDIAGUNIT, b, n, m, n, FP_ONE, A, m, B, m);

	BBLAS_LTRSM(OMPSSBLAS_LOWERTRIANG, OMPSSBLAS_TRANSP, OMPSSBLAS_NDIAGUNIT, b, n, m, n, FP_ONE, A, m, B, m);

	return 0;
}
