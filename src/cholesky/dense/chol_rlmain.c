#include "chol_rlmain.h"
#include "task_gemm.h"
#include "tasks_syrk.h"
#include "tasks_trsm.h"
#include "tasks_potrf.h"
#include "tasks_trsm.h"


#include <sys/time.h>
#include <stddef.h>
#include <stdio.h>

#include "fptype.h"


#define PRIOR_SYRK			9				
#define PRIOR_POTRF			9999
#define PRIOR_GEMM			9
#define PRIOR_TRSM			9999


int CHOL_RL(int m, int b, int t, fp_t *A, int lda) 
{

	int k;
	for ( k = 0; k < m; k += b ) {
		int doffs = k * (lda + 1);
		int bs = m-k > b ? b : m-k;

		TASK_POTRF(OMPSSLAPACK_LOWERTRIANG, bs, &A[doffs], lda, PRIOR_POTRF);

		int mk = m - k;
		int i;

		for ( i = b; i < mk; i += b ) {
			int ibs = m-i > b ? b : m-i;
			TASK_TRSM(OMPSSBLAS_RIGHT, OMPSSBLAS_LOWERTRIANG, OMPSSBLAS_TRANSP, OMPSSBLAS_NDIAGUNIT, ibs, bs, FP_ONE, &A[doffs], lda, &A[doffs + i], lda, PRIOR_TRSM);
		}

		for ( i = b; i < mk; i += b ) {
			int ibs = m-k-i > b? b : m-k-i;
			int j;
			for ( j = b; j < i; j += b ) {
				int jbs = m-j > b ? b : m-j;
				TASK_GEMM(OMPSSBLAS_NTRANSP, OMPSSBLAS_TRANSP, ibs, jbs, bs, FP_MONE, &A[doffs + i], lda, &A[doffs + j], lda, FP_ONE, &A[(k+j)*m+k+i], lda, PRIOR_GEMM);
			}

			TASK_SYRK(OMPSSBLAS_LOWERTRIANG, OMPSSBLAS_NTRANSP, ibs, bs, FP_MONE, &A[doffs + i], lda, FP_ONE, &A[(k+i)*m+k+i], lda, PRIOR_SYRK);
		}
	}

	return 0;
}


#undef PRIOR_SYRK			
#undef PRIOR_POTRF		
#undef PRIOR_GEMM	
#undef PRIOR_TRSM
