#include "lurecurs_main.h"
#include "lu_kernels.h"


int *IPIVoffs;
int IPIVi=0;
fp_t *Astart;


void lu_recurs_l(int ldim, fp_t *A, int *IPIV) {
	int skip = IPIVoffs[0];

	int j1;
	for ( j1= 1; j1 < IPIVi ; j1++ ) {
		int *ipiv = &IPIV[skip];
		int updw = IPIVoffs[j1];

		fp_t *Al = A;
		int j2;
		for ( j2=0; j2 < j1; j2++ ) {
			int aw = IPIVoffs[j2]; 
			
			task__dlaswp(skip,ldim,aw,updw,Al,ipiv);
			
			Al += aw*ldim;;
		}

		skip += IPIVoffs[j1];
	}
}


void lu_recurs_rupd(int skip, int updw, int ldim, int n, int t, fp_t *A, int *IPIV, fp_t *B) {
	if ( n <= t) {
		task__dtsrm_dgemm(skip, ldim, n, updw, n, A, IPIV, B);
	} else {
		int n2 = n>>1;
		int n1 = n - n2;

		lu_recurs_rupd(skip, updw, ldim, n1, t, A, IPIV, B);

		fp_t *Br = &B[n1 * ldim];

		lu_recurs_rupd(skip, updw, ldim, n2, t, A, IPIV, Br);
	}
}

// Submatrices of A and IPIV are passed as pointers to the first element in the first column. Skip indicates 
// the offset at which the computation starts. 
// The leading dimension (or the number of rows, same thing) plus the skip completely define the dimensions 
// of the submatrix.
void lu_recurs(int depth,int skip, int ldim, int n, int t, fp_t *A, int *IPIV) {
	if ( n <= t) {
		IPIVoffs[IPIVi++] = n;
		task__dgetf2(skip, ldim, n, A, IPIV);
	} else {
		int n2 = n>>1;
		int n1 = n - n2;

		int start1 = IPIVi;
		lu_recurs(depth+1, skip, ldim, n1, t, A, IPIV);
		int stop1 = IPIVi;

		int jj;
		fp_t *Br = &A[n1 * ldim];
		int acc=0;
		int upd;
		for ( upd=start1; upd < IPIVi; upd++ ) {
			int updw = IPIVoffs[upd]; 
			int *ipiv = &IPIV[acc];
			fp_t *Ar = &A[acc * ldim];

			lu_recurs_rupd(skip+acc, updw, ldim, n2, t, Ar, ipiv, Br); 

			acc += updw;
		}

		skip += n1;
		IPIV += n1;
		
		fp_t *Ar = &A[n1 * ldim];
		int start2 = IPIVi;
		lu_recurs(depth+1, skip, ldim, n2, t, Ar, IPIV);
	}
}
