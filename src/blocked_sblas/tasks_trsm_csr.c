#include "tasks_trsm_csr.h"


#include "blas.h"
#include "fpblas.h"
#include "fpsblas.h"
#include "tasks_potrf_csr.h"


#ifdef DOUBLE_PRECISION

#define __t_trsm_csr		task_dtrsm_csr

#else

#define __t_trsm_csr		task_strsm_csr

#endif


void __t_trsm_csr(char *trans, fp_t alpha, char *matdescra, hbmat_t *A, fp_t *X, fp_t *Y) 
{
	if ( !A->FACT ) {
		TASK_POTRF_CSR(A);
		A->FACT = 1;
	}

	int m = A->m; 
	int n = A->n;
	int *vptr = A->vptr; 
	int *vpos = A->vpos; 
	fp_t *vval = A->vval;

	SBLAS_csrsv(trans, &m, &alpha, matdescra, vval, vpos, vptr, vptr+1, X, Y);
}
