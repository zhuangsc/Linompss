#include "chol_llmain.h"

#include <stddef.h>

#include "task_gemm.h"
#include "tasks_syrk.h"
#include "tasks_trsm.h"
#include "tasks_potrf.h"
#include "tasks_trsm.h"

#include "fptype.h"
#include "fpblas.h"
#include "fpmatr.h"
#include "fplapack.h"


#define PRIOR_SYRK			9999				
#define PRIOR_POTRF			9999
#define PRIOR_GEMM			9
#define PRIOR_TRSM			9


int CHOL_LL(int m, int b, int t, fp_t *A, int lda) 
{
	int bm = b * lda;

	int coffs = 0;
	int k;
	for ( k = 0; k < m; k += b ) {
		int doffs = coffs + k;
		int bs = m-k > b ? b: m-k;

		int j;
		for ( j = 0; j < coffs; j += bm ) {
			TASK_SYRK(OMPSSBLAS_LOWERTRIANG, OMPSSBLAS_NTRANSP, bs, b, FP_MONE, &A[ j + k], lda, FP_ONE, &A[ doffs ], lda, PRIOR_SYRK);
		}

		TASK_POTRF(OMPSSLAPACK_LOWERTRIANG, bs, &A[doffs], lda, PRIOR_POTRF);

		int i;
		for ( i = k + b; i < m; i += b ) {
			int bdoffs = coffs + i;
			int ibs = m-i > b ? b : m-i;

			int j;
			for ( j=0; j<coffs; j += bm ) {
				TASK_GEMM(OMPSSBLAS_NTRANSP, OMPSSBLAS_TRANSP, ibs, b, b, FP_MONE, &A[j + i], lda, &A[j + k], lda, FP_ONE, &A[ bdoffs ], lda, PRIOR_GEMM);
			}

			TASK_TRSM(OMPSSBLAS_RIGHT, OMPSSBLAS_LOWERTRIANG, OMPSSBLAS_TRANSP, OMPSSBLAS_NDIAGUNIT, ibs, b, FP_ONE, &A[doffs], lda, &A[ bdoffs ], lda, PRIOR_TRSM);
		}
		
		coffs += bm;
	}

	return 0;
}


#undef PRIOR_SYRK			
#undef PRIOR_POTRF			
#undef PRIOR_GEMM		
#undef PRIOR_TRSM	
