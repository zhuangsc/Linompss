#include "syrk_main.h"

#include "tasks_syrk.h"
#include "task_gemm.h"

void SYRK_MAIN(ompssblas_t uplo, ompssblas_t trans, int b, int n, int k, fp_t alpha, fp_t *A, int lda, fp_t beta, fp_t *C, int ldc) 
{
	if ( !TEST_LOWER(uplo) && !TEST_UPPER(uplo) )
		fprintf(stderr, "SYRK: Parameter 1 error\n");
	if ( !TEST_NTRANSP(trans) && !TEST_TRANSP(trans) )
		fprintf(stderr, "SYRK: Parameter 2 error\n");

	if ( TEST_NTRANSP(trans) ) {
		if ( TEST_LOWER(uplo) ) {
			/* NTRANS, LOWER */
			int je = 0;
			int i;
			for ( i = 0; i < n; i += b ) {
				int bl = n - i;
				int tb = bl < b ? bl : b;

				/* all cols j != i */
				int j;
				for ( j = 0; j < je; j += b ) {
					int coffs = j * ldc + i;
					fp_t *Cb = &C[coffs];
					
					/* 											   m   n   k */
					TASK_GEMM(OMPSSBLAS_NTRANSP, OMPSSBLAS_TRANSP, tb, b, b, alpha, &A[i+0*lda], lda, &A[0*lda+j], lda, beta, Cb, ldc, 1);

					int l;
					for ( l = b; l < k; l += b ) {
						int kk = k - l;
						int dd = kk < b ? kk : b;
						TASK_GEMM(OMPSSBLAS_NTRANSP, OMPSSBLAS_TRANSP, tb, b, dd, alpha, &A[i+l*lda], lda, &A[l*lda+j], lda, FP_ONE, Cb, ldc, 1);
					}	
				}

				/* col j == i */
				int coffs = i * ldc + i;
				fp_t *Cb = &C[coffs];
				TASK_SYRK(OMPSSBLAS_LOWERTRIANG, OMPSSBLAS_NTRANSP, tb, b, alpha, &A[i], lda, beta, Cb, ldc, 1);

				int l;
				for ( l = b; l < k; l += b ) {
					int ll = k - l;
					int td = ll < b ? ll : b;
					TASK_SYRK(OMPSSBLAS_LOWERTRIANG, OMPSSBLAS_NTRANSP, tb, td, alpha, &A[i+lda*l], lda, FP_ONE, Cb, ldc, 1);
				}

				je += b;
			}
		} else {
			/* NTRANS, UPPER */
			int ie = 0;
			int j;
			for ( j = 0; j < n; j += b ) {
				int bl = n - j;
				int tb = bl < b ? bl : b;

				/* all cols j != i */
				int i;
				for ( i = 0; i < ie; i += b ) {
					int coffs = j * ldc + i;
					fp_t *Cb = &C[coffs];
					
					/* 											   m   n   k */
					TASK_GEMM(OMPSSBLAS_NTRANSP, OMPSSBLAS_TRANSP, b, tb, b, alpha, &A[i+0*lda], lda, &A[0*lda+j], lda, beta, Cb, ldc, 1);

					int l;
					for ( l = b; l < k; l += b ) {
						int kk = k - l;
						int dd = kk < b ? kk : b;
						TASK_GEMM(OMPSSBLAS_NTRANSP, OMPSSBLAS_TRANSP, b, tb, dd, alpha, &A[i+l*lda], lda, &A[l*lda+j], lda, FP_ONE, Cb, ldc, 1);
					}	
				}

				/* col j == i */
				int coffs = i * ldc + i;
				fp_t *Cb = &C[coffs];
				TASK_SYRK(OMPSSBLAS_UPPERTRIANG, OMPSSBLAS_NTRANSP, tb, b, alpha, &A[i], lda, beta, Cb, ldc, 1);

				int l;
				for ( l = b; l < k; l += b ) {
					int ll = k - l;
					int td = ll < b ? ll : b;
					TASK_SYRK(OMPSSBLAS_UPPERTRIANG, OMPSSBLAS_NTRANSP, tb, td, alpha, &A[i+lda*l], lda, FP_ONE, Cb, ldc, 1);
				}

				ie += b;
			}
		}
	} else {
		if ( TEST_LOWER(uplo) ) {
			/* TRANS, LOWER */
			int je = 0;
			int i;
			for ( i = 0; i < n; i += b ) {
				int bl = n - i;
				int tb = bl < b ? bl : b;

				/* all cols j != i */
				int j;
				for ( j = 0; j < je; j += b ) {
					int coffs = j * ldc + i;
					fp_t *Cb = &C[coffs];
					
					/* 											   m   n   k */
					TASK_GEMM(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, tb, b, b, alpha, &A[i*lda], lda, &A[j*lda], lda, beta, Cb, ldc, 1);

					int l;
					for ( l = b; l < k; l += b ) {
						int kk = k - l;
						int dd = kk < b ? kk : b;
						TASK_GEMM(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, tb, b, dd, alpha, &A[l+i*lda], lda, &A[j*lda+l], lda, FP_ONE, Cb, ldc, 1);
					}	
				}

				/* col j == i */
				int coffs = i * ldc + i;
				fp_t *Cb = &C[coffs];
				TASK_SYRK(OMPSSBLAS_LOWERTRIANG, OMPSSBLAS_TRANSP, tb, b, alpha, &A[i*lda], lda, beta, Cb, ldc, 1);

				int l;
				for ( l = b; l < k; l += b ) {
					int ll = k - l;
					int td = ll < b ? ll : b;
					TASK_SYRK(OMPSSBLAS_LOWERTRIANG, OMPSSBLAS_TRANSP, tb, td, alpha, &A[l+lda*i], lda, FP_ONE, Cb, ldc, 1);
				}

				je += b;
			}
		} else {
			/* TRANS, UPPER */
			int ie = 0;
			int j;
			for ( j = 0; j < n; j += b ) {
				int bl = n - j;
				int tb = bl < b ? bl : b;

				/* all cols j != i */
				int i;
				for ( i = 0; i < ie; i += b ) {
					int coffs = j * ldc + i;
					fp_t *Cb = &C[coffs];
					
					/* 											   m   n   k */
					TASK_GEMM(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, b, tb, b, alpha, &A[i*lda], lda, &A[j*lda], lda, beta, Cb, ldc, 1);

					int l;
					for ( l = b; l < k; l += b ) {
						int kk = k - l;
						int dd = kk < b ? kk : b;
						TASK_GEMM(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, b, tb, dd, alpha, &A[l+i*lda], lda, &A[j*lda+l], lda, FP_ONE, Cb, ldc, 1);
					}	
				}

				/* col j == i */
				int coffs = i * ldc + i;
				fp_t *Cb = &C[coffs];
				TASK_SYRK(OMPSSBLAS_UPPERTRIANG, OMPSSBLAS_TRANSP, tb, b, alpha, &A[i*lda], lda, beta, Cb, ldc, 1);

				int l;
				for ( l = b; l < k; l += b ) {
					int ll = k - l;
					int td = ll < b ? ll : b;
					TASK_SYRK(OMPSSBLAS_UPPERTRIANG, OMPSSBLAS_TRANSP, tb, td, alpha, &A[l+lda*i], lda, FP_ONE, Cb, ldc, 1);
				}

				ie += b;
			}
		}
	}
}
