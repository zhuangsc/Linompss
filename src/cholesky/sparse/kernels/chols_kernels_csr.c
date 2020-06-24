#include "chols_kernels.h"
#include "hb.h"


#ifdef SINGLE_PRECISION

#define __potrf_sparse_csr			spotrf_sparse_csr
#define __syrk_sparse_csr			ssyrk_sparse_csr
#define __gemm_sparse_csr			sgemm_sparse_csr
#define __trsm_sparse_csr			strsm_sparse_csr

#else

#define __potrf_sparse_csr			dpotrf_sparse_csr
#define __syrk_sparse_csr			dsyrk_sparse_csr
#define __gemm_sparse_csr			dgemm_sparse_csr
#define __trsm_sparse_csr			dtrsm_sparse_csr

#endif

/* A in CSR, lower-triangular */
void __potrf_sparse_csr(hbmat_t* A) {

 	/* Check if the input matrix has properly set */
	if ( A->vptr == NULL ) {
		hyper_sym_csr_task2(A);
	}

	int m = A->m;
	int* vptr = A->vptr; int* vpos = A->vpos; 
	fp_t* vval = A->vval;

	int d;
	for ( d=0 ; d<m; ++d ) {
		int dstart = vptr[d];
		int dend = vptr[d+1];

		int diagidx = dend - 1;
		fp_t diagel = vval[diagidx];
		/* compute element on diagonal (d,d) */
		int j;
		for ( j = dstart; j < diagidx; ++j ) {
			diagel -= vval[j]*vval[j];
		}
		diagel = sqrt(diagel);
		vval[diagidx] = diagel;


		int i;
		for ( i = d+1; i < m; ++i ) {
			int jx = vptr[i];
			int j = vpos[jx];
			int jend = vptr[i+1];

			int djx = dstart;
			int dj = vpos[djx];

			fp_t acc = 0;
			while ( dj < d && djx < dend && j < d && jx < jend ) {
				int djadv = dj < j; 
				int jadv  = j  < dj; 

				jx  += jadv? 1 : 0;
				djx += djadv? 1 : 0;
				dj = vpos[djx];
			 	j  = vpos[ jx];

				if ( dj == j && j < d && djx < dend && jx < jend ) {
					acc += vval[djx++] * vval[jx++];
					dj = vpos[djx];
			 		j  = vpos[ jx];
				}
			}

			while ( vpos[jx] < d && jx < jend ) { ++jx; }

			if ( jx < jend && vpos[jx]==d ) {
				fp_t subdval = vval[jx] - acc;
				vval[jx] = subdval / diagel;
			}
		}
	}
}


/* A and C in CSR, lower-triangular */
void __syrk_sparse_csr(hbmat_t* A, hbmat_t* C) {
#if USE_MKL
	 /* Check if the input matrix has properly set */
	if ( A->vptr == NULL )
		hyper_sym_csr_task2(A);
	if ( C->vptr == NULL )
		hyper_sym_csr_task2(C);

//	hb_sanity_check("A", A, 0);
//	hb_sanity_check("C", C, 0);
	int n = A->n; int m = A->m;
	int* vptr = A->vptr; int* vpos = A->vpos; fp_t* vval = A->vval;
	int* vptr_c = C->vptr; int* vpos_c = C->vpos; fp_t* vval_c = C->vval;

	fp_t* peela = malloc(2 * m * sizeof(fp_t));
	fp_t* peelc = &peela[m];

	int i;
	for ( i = 0; i < m; ++i ) {
		array_s2d(A, peela, i);

#ifdef SINGLE_PRECISION
		mkl_cspblas_scsrgemv("N", &m, vval, vptr, vpos, peela, peelc);
#else
		mkl_cspblas_dcsrgemv("N", &m, vval, vptr, vpos, peela, peelc);
#endif

		int k;
		for ( k = vptr_c[i]; k < vptr_c[i+1]; ++k ) {
			int col_pos = vpos_c[k];
			vval_c[k] -= peelc[col_pos];
		}

		array_clear2(A, peela, i);
	}

	free(peela); 
#endif
}


/* A, B and C in CSR */
/* C = C - A * B^T, or */
/* C^T = C^T - B * A^T */
void __gemm_sparse_csr(hbmat_t* A, hbmat_t* B, hbmat_t* C) {
#if USE_MKL
 	/* Check if the input matrix has properly set */
	if ( A->vptr == NULL ) {
		hyper_sym_csr_task2(A);
	}

	if ( C->vptr == NULL ) {
		hyper_sym_csr_task2(C);
	}

	int m = B->m; int n = B->n;
	int* vptr = B->vptr; int* vpos = B->vpos; fp_t* vval = B->vval;
	int* vptr_c = C->vptr; int* vpos_c = C->vpos; fp_t* vval_c = C->vval;

	fp_t* peelb = malloc(2*m*sizeof(fp_t));
	fp_t* peelc = &peelb[m];

	int i;
	for ( i = 0; i < m; ++i ) {
		array_clear(peelb, m);
		array_s2d(A, peelb, i);

#ifdef SINGLE_PRECISION
		mkl_cspblas_scsrgemv("N", &m, vval, vptr, vpos, peelb, peelc);
#else
		mkl_cspblas_dcsrgemv("N", &m, vval, vptr, vpos, peelb, peelc);
#endif

		int k;
		for ( k = vptr_c[i]; k < vptr_c[i+1]; k++ ) {
			int col_pos = vpos_c[k];
			vval_c[k] -= peelc[col_pos];
		}
	}

	free(peelb); 
#endif
}


/* A CSR, C CSR */
/* C = C \ A */
void __trsm_sparse_csr(hbmat_t* A, hbmat_t* C){
#if USE_MKL
	/* Check if the input matrix has properly set */
	if ( C->vptr == NULL ) {
		hyper_sym_csr_task2(C);
	}

	int n = A->n; int m = A->m;
	int* vptr = A->vptr; int* vpos = A->vpos; fp_t* vval = A->vval;

	fp_t* peelb = calloc( 2 * m, sizeof(fp_t));
	fp_t* peelc = &peelb[m];
	fp_t alpha = 1;

	int i;
	for ( i = 0; i < n; ++i ) {
		array_s2d(C, peelb, i);

#ifdef SINGLE_PRECISION
		mkl_scsrsv("N", &n, &alpha, "TLNC", vval, vpos, vptr, vptr+1, peelb, peelc);
#else
		mkl_dcsrsv("N", &n, &alpha, "TLNC", vval, vpos, vptr, vptr+1, peelb, peelc);
#endif

		array_d2s(C, peelc, i);
		array_clear2(C, peelb, i);
	}

	free(peelb); 
#endif
}
