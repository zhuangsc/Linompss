#include "chol_rlmain.h"
#include "task_gemm.h"
#include "tasks_syrk.h"
#include "tasks_trsm.h"
#include "tasks_potrf.h"
#include "tasks_trsm.h"


#include <stddef.h>

#include "fptype.h"

#include "fpblas.h"


int CHOL_HRL(int mt, int b, int t, fp_t **Ah) {
	int k;
	for ( k = 0; k < mt ;k++ ) {

		TASK_POTRF(OMPSSLAPACK_LOWERTRIANG, b, Ah[k * mt + k], b, 1);

		int i;
		for (i = k+1; i < mt; i++) {
			TASK_TRSM(OMPSSBLAS_RIGHT, OMPSSBLAS_LOWERTRIANG, OMPSSBLAS_TRANSP, OMPSSBLAS_NDIAGUNIT, b, b, FP_ONE, Ah[k * mt + k], b, Ah[k * mt + i], b, 1);
		}

		for (i = k+1; i < mt; i++) {
			int j;
			for (j=k+1; j<i; j++) {
				TASK_GEMM(OMPSSBLAS_NTRANSP, OMPSSBLAS_TRANSP, b, b, b, FP_MONE, Ah[k * mt + i], b, Ah[k * mt + j], b, FP_ONE, Ah[ j * mt + i], b, 1);
			}
			TASK_SYRK(OMPSSBLAS_LOWERTRIANG, OMPSSBLAS_NTRANSP, b, b, FP_MONE, Ah[k * mt + i], b, FP_ONE, Ah[i * mt + i], b, 1);
		}
	}

	return 0;
}
