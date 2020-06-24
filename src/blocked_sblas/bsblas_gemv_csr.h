#ifndef __BSBLAS_GEMV_CSR_H__
#define __BSBLAS_GEMV_CSR_H__


#include "fptype.h"
#include "hb.h"
#include "async.h"
#include "async_struct.h"
#include "tasks_gemv_csr.h"
#include "bblas_copy.h"


#ifdef SINGLE_PRECISION

#define BSBLAS_NODIAG_CP_GEMV_CSR			bsblas_nd_scpgemv_csr
#define BSBLAS_NODIAG_CP_GEMV_CSR_PRED		bsblas_nd_scpgemv_csr_pred
#define BSBLAS_CP_GEMV_CSR					bsblas_scpgemv_csr

#else

#define BSBLAS_NODIAG_CP_GEMV_CSR			bsblas_nd_dcpgemv_csr
#define BSBLAS_NODIAG_CP_GEMV_CSR_PRED		bsblas_nd_dcpgemv_csr_pred
#define BSBLAS_CP_GEMV_CSR					bsblas_dcpgemv_csr

#endif

/*
 * y = b - A * x 
 * skip the diagonal blocks
 */
static inline void __attribute__((always_inline)) bsblas_nd_scpgemv_csr(hbmat_t *A, int bs, \
		int alpha, float *X, float beta, float *B, float *Y) {
	int N = A->n;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	hbmat_t **vval = A->vval;
	int I;
	BBLAS_COPY(OMPSS_PRIOR_DFLT, bs, 1, A->orig->m, 1, B, Y);
	for ( I = 0; I < N; ++I ) {
		int J = vptr[I];
		for ( ; J < vptr[I+1]; ++J ) {
			if ( vpos[J] != I ) {
				TASK_GEMV_CSR(vval[J], "N", alpha, "GLNC", beta, &(X[vpos[J] * bs]), &(Y[I * bs]));
			}
		}
	}
}

static inline void __attribute__((always_inline)) bsblas_nd_dcpgemv_csr(hbmat_t *A, int bs, \
		double alpha, double *X, double beta, double *B, double *Y) {
	int N = A->n;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	hbmat_t **vval = A->vval;
	BBLAS_COPY(OMPSS_PRIOR_DFLT, bs, 1, A->orig->m, 1, B, Y);
	int I;
	for ( I = 0; I < N; ++I ) {
		int J = vptr[I];
		for ( ; J < vptr[I+1]; ++J ) {
			if ( vpos[J] != I ) {
				TASK_GEMV_CSR(vval[J], "N", alpha, "GLNC", beta, &(X[vpos[J] * bs]), &(Y[I * bs]));
			}
		}
	}
}


static inline async_stat_t __attribute__((always_inline)) bsblas_nd_scpgemv_csr_pred(int it, hbmat_t *A, int bs, \
		float alpha, float *X, float beta, float *B, float *Y, float threshold, float *cdot, async_t *dotsync, fp_t div) {
	int N = A->n;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	hbmat_t **vval = A->vval;
	async_stat_t status = STAT_AHEAD;
	BBLAS_COPY(OMPSS_PRIOR_DFLT, bs, 1, A->orig->m, 1, B, Y);
	int I;
	for ( I = 0; I < N; ++I ) {
		status = ASYNC_CONV(it, threshold, cdot, div, dotsync, 0);
		if ( status == STAT_CONVERGED ) {
			return STAT_CONVERGED;
		}
		int J = vptr[I];
		for ( ; J < vptr[I+1]; ++J ) {
			if ( vpos[J] != I ) {
				TASK_GEMV_CSR(vval[J], "N", alpha, "GLNC", beta, &(X[vpos[J] * bs]), &(Y[I * bs]));
			}
		}
	}
	return status;
}

static inline async_stat_t __attribute__((always_inline)) bsblas_nd_dcpgemv_csr_pred(int it, hbmat_t *A, int bs, \
		double alpha, double *X, double beta, double *B, double *Y, double threshold, double *cdot, async_t *dotsync, fp_t div) {
	int N = A->n;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	hbmat_t **vval = A->vval;
	async_stat_t status;
	BBLAS_COPY(OMPSS_PRIOR_DFLT, bs, 1, A->orig->m, 1, B, Y);
	int I;
	for ( I = 0; I < N; ++I ) {
		status = ASYNC_CONV(it, threshold, cdot, div, dotsync, 0);
		if ( status == STAT_CONVERGED ) {
			return STAT_CONVERGED;
		}
		int J = vptr[I];
		for ( ; J < vptr[I+1]; ++J ) {
			if ( vpos[J] != I ) {
				TASK_GEMV_CSR(vval[J], "N", alpha, "GLNC", beta, &(X[vpos[J] * bs]), &(Y[I * bs]));
			}
		}
	}
	return status;
}


/*
 * y = b - A * x 
 */

static inline void __attribute__((always_inline)) bsblas_scpgemv_csr(hbmat_t *A, int bs, float alpha, float *X, float beta, float *B, float *Y) 
{
	int N = A->n;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	hbmat_t **vval = A->vval;
	BBLAS_COPY(OMPSS_PRIOR_DFLT, bs, 1, A->orig->m, 1, B, Y);
	int I;
	for ( I = 0; I < N; ++I ) {
		int J = vptr[I];
		TASK_GEMV_CSR(vval[J], "N", alpha, "GLNC", beta, &(X[vpos[J] * bs]), &(Y[I * bs]));
		++J;
		for ( ; J < vptr[I+1]; ++J ) {
			TASK_GEMV_CSR(vval[J], "N", alpha, "GLNC", FP_ONE, &(X[vpos[J] * bs]), &(Y[I * bs]));
		}
	}
}

static inline void __attribute__((always_inline)) bsblas_dcpgemv_csr(hbmat_t *A, int bs, double alpha, double *X, double beta, double *B, double *Y) 
{
	int N = A->n;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	hbmat_t **vval = A->vval;
	BBLAS_COPY(OMPSS_PRIOR_DFLT, bs, 1, A->orig->m, 1, B, Y);
	int I;
	for ( I = 0; I < N; ++I ) {
		int J = vptr[I];
		
		TASK_GEMV_CSR(vval[J], "N", alpha, "GLNC", beta, &(X[vpos[J] * bs]), &(Y[I * bs]));
		++J;
		for ( ; J < vptr[I+1]; ++J ) {
			TASK_GEMV_CSR(vval[J], "N", alpha, "GLNC", FP_ONE, &(X[vpos[J] * bs]), &(Y[I * bs]));
		}
	}
}

#endif
