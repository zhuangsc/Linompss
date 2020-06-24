#include "lu_main.h"
#include "lu_kernels.h"


void LU_MAIN(int m,int n,int b, fp_t *A, int lda, int *IPIV) 
{
	int j;
	for( j=0; j<n; j+=b) {
		int dimC = n-j > b ? b : n-j;
		int dimR = m - j;
		int aoffs = j * lda + j;
		int *ip = &IPIV[j];

		TASK_GETRF(j, dimR, dimC, &A[j*lda], lda, ip);
		int ips = min(dimR, dimC);
		int jj;
		for( jj=0; jj<j; jj+=b){
			TASK_LASWP(j, dimC, &A[jj*lda], lda, 1, ips, ip, 1);
		}
		int dimM = dimC;
		for( jj=j+dimM; jj<n; jj+=dimM ) {
			dimM = n-jj > b ? b : n-jj;
			TASK_TRSM_GEMM(j, &A[j*lda], &A[jj*lda], dimM, dimC, dimR, ips, ip, lda);
		}
	}
} 
