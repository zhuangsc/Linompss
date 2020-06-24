#include "lull_main.h"
#include "lu_kernels.h"


void LULL_MAIN(int m, int n, int b, fp_t *A, int lda, int *IPIV) 
{
	int j;
	for ( j = 0; j<n; j+=b ) {

		int dimC = n-j > b ? b : n-j;
		int dimR = m - j;
		int aoffs = j * lda + j;
		int ips = min(dimR, dimC);
		int ip = &IPIV[j];
		
		int jj;
		for ( jj=0; jj<j; jj+=b ) {
			ip = &IPIV[jj];
			TASK_TRSM_GEMM_LL( jj, &A[jj*lda], &A[j*lda], b, dimC, m-jj, ips, ip, lda);
		}
	
		ip = &IPIV[j];
		TASK_GETRF( j, dimR, dimC, &A[j*lda], lda, ip);

		for (jj = 0; jj < j; jj+=b) {
			TASK_LASWP( j, b, &A[jj*lda], lda, 1, ips, ip, 1);
		}
	}
}
