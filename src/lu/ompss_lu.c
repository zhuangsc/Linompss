#include "ompss_lu.h"
#include "fptype.h"
#include "fpblas.h"
#include "tasks_trsm.h"
#include "task_gemm.h"
#include "lu_main.h"
#include "lull_main.h"
#include "lu_kernels.h"

void ompss_sgetrf(int m, int n, int b, float *A, int lda, int *IPIV)
{
	slu_main( m, n, b, A, lda, IPIV);
}

void ompss_dgetrf(int m, int n, int b, double *A, int lda, int *IPIV)
{
	dlu_main( m, n, b, A, lda, IPIV);
}


void ompss_sgetrf_ll(int m, int n, int b, float *A, int lda, int *IPIV)
{
	slull_main( m, n, b, A, lda, IPIV);
}

void ompss_dgetrf_ll(int m, int n, int b, double *A, int lda, int *IPIV)
{
	dlull_main( m, n, b, A, lda, IPIV);
}

#if 0
void ompss_sgetrf2(int m, int n, int b, float *A, int lda, int *IPIV)
{
	float MONE = -1.0;
	float ONE = 1.0;
	int j;
	for( j=0; j<n; j+=b) {
		int dimC = n-j > b ? b : n-j;
		int dimR = m - j;
		int skip = j;
		int aoffs = j * lda + j;
		int *ip = &IPIV[j];

		task_sgetrf(j, dimR, dimC, &A[j*lda], lda, ip);
		int ips = min(dimR, dimC);
		int jj;
		for( jj=0; jj<j; jj+=b){
			int offs = jj * lda + j;
			task_slaswp(j, dimC, &A[jj*lda], lda, 1, ips, ip, 1);
		}
		int dimM = dimC;
		for( jj=j+dimM; jj<n; jj+=dimM ) {
			dimM = n-jj > b ? b : n-jj;
			task_strsm_sgemm(j, &A[j*lda], &A[jj*lda], dimM, dimC, dimR, ips, ip, lda);
		}
	}
}

void ompss_dgetrf2(int m, int n, int b, double *A, int lda, int *IPIV)
{
	double MONE = -1.0;
	double ONE = 1.0;
	int j;
	for( j=0; j<n; j+=b) {
		int dimC = n-j > b ? b : n-j;
		int dimR = m - j;
		int skip = j;
		int aoffs = j * lda + j;
		int *ip = &IPIV[j];

		task_dgetrf(j, dimR, dimC, &A[j*lda], lda, ip);
		int ips = min(dimR, dimC);
		int jj;
		for( jj=0; jj<j; jj+=b){
			int offs = jj * lda + j;
			task_dlaswp(j, dimC, &A[jj*lda], lda, 1, ips, ip, 1);
		}
		int dimM = dimC;
		for( jj=j+dimM; jj<n; jj+=dimM ) {
			dimM = n-jj > b ? b : n-jj;
			task_dtrsm_dgemm(j, &A[j*lda], &A[jj*lda], dimM, dimC, dimR, ips, ip, lda);
		}
	}
}
#endif
