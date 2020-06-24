#include "chols_kernels.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif 		

#include "hbconvrt.h"
#include "array.h"

#ifdef SINGLE_PRECISION

#define __potrf_sparse_csc			spotrf_sparse_csc
#define __syrk_sparse_csc			ssyrk_sparse_csc
#define __gemm_sparse_csc			sgemm_sparse_csc
#define __trsm_sparse_csc			strsm_sparse_csc

#else

#define __potrf_sparse_csc			dpotrf_sparse_csc
#define __syrk_sparse_csc			dsyrk_sparse_csc
#define __gemm_sparse_csc			dgemm_sparse_csc
#define __trsm_sparse_csc			dtrsm_sparse_csc

#endif

void __potrf_sparse_csc(hbmat_t* A) {
#if USE_MKL
	int n = A->n;
	int* vptr = A->vptr; int* vpos = A->vpos; 
	fp_t* vval = A->vval;

	int J;
	for ( J=0 ; J < n; ++J ) {
		//l22 = sqrt (a22 - l12*l12_t)
		fp_t* a22 = &(vval[vptr[J]]);
		fp_t* l22 = a22;
		int p_a22 = vpos[vptr[J]];

		int j;
		for ( j=0; j<J; j++ ) {
			int i;
			for ( i=vptr[j]; i<vptr[j+1]; i++ ) {
				if ( vpos[i] == p_a22 ) {
					*a22 -= vval[i]*vval[i];
				}
			}
		}
		*l22 = sqrt(*a22);
		// l32 = (a32-L32*l12)/l22
		int i;
		for ( i = vptr[J]+1; i < vptr[J+1]; i++ ) {
			fp_t* a32 = &(vval[i]);
			fp_t* l32 = a32;
			fp_t l31 = 0; 
			fp_t l12 = 0;
			fp_t l31l12 = 0;
			int p_a32 = vpos[i];

			int j;
			for ( j = 0; j < J; j++ ) {
				l31=l12=0;
				int k;
				for ( k=vptr[j]; k<vptr[j+1]; k++) {
					if (vpos[k] < p_a22)
						continue;
					if ( vpos[k] > p_a32 )
						break;
					if ( vpos[k] == p_a22 ) {
						l12 = vval[k];
					}
					if ( vpos[k] == p_a32 ) {
						l31 = vval[k]; 
						break;
					}
				}
				l31l12 += l31*l12;
			}
			*l32 = (*a32 - l31l12)/(*l22);
		}
	}
#endif
}

void __syrk_sparse_csc(hbmat_t* A, hbmat_t* C) {
#if USE_MKL
	int n = A->n; int m = A->m;
	int* vptr = A->vptr;
	int* vpos = A->vpos;
	fp_t* vval = A->vval;
	fp_t* peela = malloc(m*sizeof(fp_t));
	fp_t* peelc = malloc(m*sizeof(fp_t));
	char* trans = "N";
	fp_t alpha = -1; 
	fp_t beta = 1;
	char* matdescra = "GLNC";
	hbmat_t* A_csr = csc2csr(A);

	int i;
	for ( i = 0; i < n; i++ ) {
		array_clear(peela, m);
		array_clear(peelc, m);
		array_s2d(A_csr, peela, i);
		array_s2d(C, peelc, i);
		mkl_dcscmv(trans, &m, &n, &alpha, matdescra, vval, vpos, vptr, vptr+1, peela, &beta, peelc);
		array_d2s(C, peelc, i);
	}

	free(peela); free(peelc);
	hb_free(A_csr);
#endif
}

void __gemm_sparse_csc(hbmat_t* A, hbmat_t* B, hbmat_t* C) {
#if USE_MKL
	hbmat_t* B_csr = csc2csr(B);
	int m = A->m; int n = A->n;
	int* vptr = A->vptr;
	int* vpos = A->vpos;
	fp_t* vval = A->vval;
	fp_t* peelb = malloc(m*sizeof(fp_t));
	fp_t* peelc = malloc(m*sizeof(fp_t));
	char* trans = "N";
	fp_t alpha = -1;
	fp_t beta = 1;
	char* matdescra = "GLNC";

	int i;
	for ( i = 0; i < n; i++ ) {
		array_clear(peelb, m);
		array_clear(peelc, m);
		array_s2d(B_csr, peelb, i);
		array_s2d(C, peelc, i);
		mkl_dcscmv(trans, &m, &n, &alpha, matdescra, vval, vpos, vptr, vptr+1, peelb, &beta, peelc);
		array_d2s(C, peelc, i);
	}

	free(peelb); free(peelc);
	hb_free(B_csr);
#endif
}

void __trsm_sparse_csc(hbmat_t* A, hbmat_t* C) {
#if USE_MKL
	hbmat_t* C_csr = csc2csr(C);
	hbmat_t* C_hb;
	int n = A->n; int m = A->m;
	int* vptr = A->vptr;
	int* vpos = A->vpos;
	fp_t* vval = A->vval;
	fp_t* peelb = malloc(m*sizeof(fp_t));
	fp_t* peelc = malloc(m*sizeof(fp_t));
	char* trans = "N";
	fp_t alpha = 1;
	char* matdescra = "TLNC";

	int i;
	for ( i = 0; i < n; i++ ) {
		array_clear(peelc, m);
		array_clear(peelb, m);
		array_s2d(C_csr, peelb, i);
		array_s2d(C_csr, peelc, i);
		mkl_dcscsv(trans, &n, &alpha, matdescra, vval, vpos, vptr, vptr+1, peelb, peelc);
		array_d2s(C_csr, peelc, i);
	}

	int job[6] = {0,0,0,1,1,1};
	int info;
	mkl_dcsrcsc(job,&n,C_csr->vval,C_csr->vpos,C_csr->vptr,C->vval,C->vpos,C->vptr,&info);

	free(peelb); free(peelc);
	hb_free(C_csr);
#endif
}

