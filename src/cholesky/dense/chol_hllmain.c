#include "chol_llmain.h"

#include <sys/time.h>
#include <stddef.h>

#include "task_gemm.h"
#include "tasks_syrk.h"
#include "tasks_trsm.h"
#include "tasks_potrf.h"
#include "tasks_trsm.h"


#include "fptype.h"
#include "fpblas.h"


int CHOL_HLL(int mt, int b, int t, fp_t **Ah) {
	int k;
	for ( k = 0; k < mt; k++ ) {
		int j;
		for (j=0; j<k; j++) {
			TASK_SYRK(OMPSSBLAS_LOWERTRIANG, OMPSSBLAS_NTRANSP, b, b, FP_MONE, Ah[j * mt + k], b, FP_ONE, Ah[k * mt + k], b, 1);
		}

		TASK_POTRF(OMPSSLAPACK_LOWERTRIANG, b, Ah[k * mt + k], b, 1);

		int i;
		for (i = k+1; i < mt; i++) {
			int j;
			for (j=0; j<k; j++) {
				TASK_GEMM(OMPSSBLAS_NTRANSP, OMPSSBLAS_TRANSP, b, b, b, FP_MONE, Ah[j * mt + i], b, Ah[j * mt + k], b, FP_ONE, Ah[ k * mt + i], b, 1);
			}

			TASK_TRSM(OMPSSBLAS_RIGHT, OMPSSBLAS_LOWERTRIANG, OMPSSBLAS_TRANSP, OMPSSBLAS_NDIAGUNIT, b, b, FP_ONE, Ah[k * mt + k], b, Ah[k * mt + i], b, 1);
		}

#if 0
		for (i = k+1; i < mt; i++) {
		      TRSM_TASK(b, t, Ah[k * mt + k], Ah[k * mt + i] );
		}
#endif
	}

	return 0;
}
