#include "tasks_gemv_csr.h"


#ifdef DOUBLE_PRECISION
#define __t_gemv_csr		task_dgemv_csr
#define __t_cpgemv_csr		task_dcpgemv_csr
#else
#define __t_gemv_csr		task_sgemv_csr
#define __t_cpgemv_csr		task_scpgemv_csr
#endif

void __t_gemv_csr(hbmat_t *A, char *trans, fp_t alpha, char *matdescra, fp_t beta, fp_t *X, fp_t *Y)
{
	int m = A->m; int n = A->n;
	int *vptr = A->vptr; int *vpos = A->vpos; fp_t *vval = A->vval;

	SBLAS_csrmv(trans, m, n, alpha, matdescra, vval, vpos, vptr, vptr+1, X, beta, Y);
}

void __t_cpgemv_csr(hbmat_t *A, char *trans, fp_t alpha, char *matdescra, fp_t beta, fp_t *X, fp_t *B, fp_t *Y) {
	int m = A->m; int n = A->n;
	int *vptr = A->vptr; int *vpos = A->vpos; fp_t *vval = A->vval;

	int i;
	for ( i = 0; i < m; ++i ) {
		Y[i] = B[i];
	}

	SBLAS_csrmv(trans, m, n, alpha, matdescra, vval, vpos, vptr, vptr+1, X, beta, Y);
}
