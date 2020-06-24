#include "lu_kernels.h"

#include <stdlib.h>

void TASK_GETRF(int skip, int m, int n, fp_t *A, int lda, int *IPIV)
{
	A += skip;
	int info;
	LAPACK_getrf(m, n, A, lda, IPIV, &info);
}

void TASK_LASWP(int skip, int n, fp_t *A, int lda, int k1, int k2, int *IPIV, int inc)
{
	A += skip;
	LAPACK_laswp(n, A, lda, k1, k2, IPIV, inc);
}

void TASK_TRSM_GEMM(int skip, fp_t *A, fp_t *B, int dimM, int dimC, int dimR, int ips, int *ip, int lda)
{

#ifdef SINGLE_PRECISION
	slaswp_(&dimM, &B[skip], &lda, &I_ONE, &ips, ip, &I_ONE);
#endif
#ifdef DOUBLE_PRECISION
	dlaswp_(&dimM, &B[skip], &lda, &I_ONE, &ips, ip, &I_ONE);
#endif
//	LAPACK_laswp(dimM, &B[skip], lda, IONE, ips, ip, IONE);
	BLAS_trsm(OMPSSBLAS_LEFT, OMPSSBLAS_LOWERTRIANG, OMPSSBLAS_NTRANSP, OMPSSBLAS_DIAGUNIT, dimC, dimM, FP_ONE, &A[skip], lda, &B[skip], lda);
	int dimdiff = dimR - dimC;
	BLAS_gemm(OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, dimdiff, dimM ,dimC, FP_MONE, &A[skip+dimC], lda, &B[skip], lda, FP_ONE, &B[skip+dimC], lda);

}

void TASK_TRSM_GEMM_LL(int skip, fp_t *A, fp_t *B, int dimM, int dimC, int dimR, int ips, int *ip, int lda)
{
	LAPACK_laswp(dimM, &B[skip], lda, I_ONE, ips, ip, I_ONE);
	BLAS_trsm(OMPSSBLAS_LEFT, OMPSSBLAS_LOWERTRIANG, OMPSSBLAS_NTRANSP, OMPSSBLAS_DIAGUNIT, dimM, dimC, FP_ONE, &A[skip], lda, &B[skip], lda);
	int dimdiff = dimR - dimM;
	BLAS_gemm(OMPSSBLAS_NTRANSP, OMPSSBLAS_NTRANSP, dimdiff, dimC ,dimM, FP_MONE, &A[skip+dimM], lda, &B[skip], lda, FP_ONE, &B[skip+dimM], lda);
}
