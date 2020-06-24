#ifndef __BSBLAS_GEMM_CSR_H__
#define __BSBLAS_GEMM_CSR_H__


#include "fptype.h"
#include "hb.h"
#include "async.h"


#define LIBSBBLAS_EXPORT __attribute__((__visibility__("default")))


#ifdef SINGLE_PRECISION
#define BSBLAS_CP_GEMM_CSR				bsblas_scpgemv_csr
#else
#define BSBLAS_CP_GEMM_CSR				bsblas_dcpgemv_csr
#endif


static inline void __attribute__((always_inline)) bsblas_scpgemv_csr(hbmat_t *A, int bs, float alpha, float *X, float beta, float *B, float *Y) 
{
	int N = A->n;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	hbmat_t **vval = A->vval;
	int I;
	for ( I = 0; I < N; ++I ) {
		int J = vptr[I];
		TASK_CPGEMM_CSR(vval[J], "N", alpha, "GLNC", beta, &(X[vpos[J] * bs]), &(B[I * bs]), &(Y[I * bs]));
		++J;
		for ( ; J < vptr[I+1]; ++J ) {
			TASK_GEMM_CSR(vval[J], "N", alpha, "GLNC", beta, &(X[vpos[J] * bs]), &(Y[I * bs]));
		}
	}
}

static inline void __attribute__((always_inline)) bsblas_dcpgemv_csr(hbmat_t *A, int bs, double alpha, double *X, double beta, double *B, double *Y) 
{
	int N = A->n;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	hbmat_t **vval = A->vval;
	int I;
	for ( I = 0; I < N; ++I ) {
		int J = vptr[I];
		TASK_CPGEMM_CSR(vval[J], "N", alpha, "GLNC", beta, &(X[vpos[J] * bs]), &(B[I * bs]), &(Y[I * bs]));
		++J;
		for ( ; J < vptr[I+1]; ++J ) {
			TASK_GEMM_CSR(vval[J], "N", alpha, "GLNC", beta, &(X[vpos[J] * bs]), &(Y[I * bs]));
		}
	}
}

#endif
