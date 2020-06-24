#ifndef __BBLAS_TRSM_H__
#define __BBLAS_TRSM_H__


#include "fptype.h"
#include "fpmatr.h"
#include "tasks_trsm.h"
#include "task_gemm.h"

#include "matfprint.h"


#ifdef SINGLE_PRECISION

#define BBLAS_LTRSM bblas_sltrsm

#endif

#ifdef DOUBLE_PRECISION

#define BBLAS_LTRSM bblas_dltrsm

#endif


/* Solve A * X = B, with A triangular */
/* A divided in blocks of bm x bm, B in blocks of bm x bn */
static inline void __attribute__((always_inline)) bblas_sltrsm(const ompssblas_t uplo, const ompssblas_t trans, const ompssblas_t diag, 
		int bm, int bn, int m, int n, float alpha, float *A, int lda, float *B, int ldb) {

	if ( TEST_LOWER(uplo) ) {
		if ( TEST_NTRANSP(trans) ) {
			/* lower, ntrans */
			task_strsm(OMPSSBLAS_LEFT, uplo, trans, diag, bm, bn, alpha, &A[0], lda, &B[0], ldb, 1); 

			int i;
			for ( i=bm; i<m; i+=bm ) {
				int dimA = m-i > bm ? bm : m-i;
				task_sgemm(trans, OMPSSBLAS_NTRANSP, dimA, bn, bm, FP_MONE, &A[i], lda, &B[0], ldb, alpha, &B[i], ldb, 1);

				int offs = i + bm * lda; 
				int j;
				for ( j=bm; j<i; j+=bm ) {
					int dimB = m-j > bm ? bm : m-j;
					task_sgemm(trans, OMPSSBLAS_NTRANSP, dimA, bn, dimB, FP_MONE, &A[offs], lda, &B[j], ldb, FP_ONE, &B[i], ldb, 1);
					offs += bm * lda;
				}

				task_strsm(OMPSSBLAS_LEFT, uplo, trans, diag, dimA, bn, FP_ONE, &A[offs], lda, &B[i], ldb, 1); 
			}
		} else {
			/* lower, trans */
			int bR = m%bm > 0 ? m%bm : bm;
			int boffs = m - bR;
			int aoffs = boffs * lda + boffs;
			int step = bm;

			task_strsm(OMPSSBLAS_LEFT, uplo, trans, diag, bR, bn, alpha, &A[aoffs], lda, &B[boffs], ldb, 1); 

			int j;
			for ( j = boffs-bm; j >= 0; j-=bm, step+=bm ) {
				aoffs -= bm * lda;
				task_sgemm(trans, OMPSSBLAS_NTRANSP, bm, bn, bR, FP_MONE, &A[aoffs], lda, &B[boffs], ldb, alpha, &B[j], ldb, 1);

				int i;
				for ( i = bm; i < step; i+=bm ) {
					task_sgemm(trans, OMPSSBLAS_NTRANSP, bm, bn, bm, FP_MONE, &A[aoffs-i], lda, &B[boffs-i], ldb, FP_ONE, &B[j], ldb, 1);
				}

				task_strsm(OMPSSBLAS_LEFT, uplo, trans, diag, bm, bn, FP_ONE, &A[j*lda+j], lda, &B[j], ldb, 1); 
			}
		}
	} else {
		if ( TEST_NTRANSP(trans) ) {
			/* upper, ntrans */
			int bR = m%bm > 0 ? m%bm : bm;
			int boffs = m - bR;
			int aoffs = boffs * lda + boffs;
			int step = bm;

			task_strsm(OMPSSBLAS_LEFT, uplo, trans, diag, bR, bn, alpha, &A[aoffs], lda, &B[boffs], ldb, 1);

			int i;
			for ( i = boffs-bm; i >= 0; i-=bm, step+=bm ) {
				aoffs -= bm;
				task_sgemm(trans, OMPSSBLAS_NTRANSP, bm, bn, bR, FP_MONE, &A[aoffs], lda, &B[boffs], ldb, alpha, &B[i], ldb, 1);

				int offs = aoffs;
				int j;
				for ( j = bm; j < step; j+=bm ) {
					offs -= bm * lda;
					task_sgemm(trans, OMPSSBLAS_NTRANSP, bm, bn, bm, FP_MONE, &A[offs], lda, &B[boffs-j], ldb, FP_ONE, &B[i], ldb, 1);
				}
				task_strsm(OMPSSBLAS_LEFT, uplo, trans, diag, bm, bn, FP_ONE, &A[i*lda+i], lda, &B[i], ldb, 1);
			}
		} else {
			/* upper, trans */
			task_strsm(OMPSSBLAS_LEFT, uplo, trans, diag, bm, bn, alpha, &A[0], lda, &B[0], ldb, 1); 

			int j;
			for ( j=bm; j<m; j+=bm ) {
				int dimA = m-j > bm ? bm : m-j;
				task_sgemm(trans, OMPSSBLAS_NTRANSP, dimA, bn, bm, FP_MONE, &A[j*lda], lda, &B[0], ldb, alpha, &B[j], ldb, 1);

				int i;
				for ( i=bm; i<j; i+=bm ) {
					int dimB = m-i > bm ? bm : m-i;
					task_sgemm(trans, OMPSSBLAS_NTRANSP, dimA, bn, dimB, FP_MONE, &A[j*lda+i], lda, &B[i], ldb, FP_ONE, &B[j], ldb, 1);
				}
				task_strsm(OMPSSBLAS_LEFT, uplo, trans, diag, dimA, bn, FP_ONE, &A[j*lda+i], lda, &B[j], ldb, 1); 
			}
		}
	}
}

static inline void __attribute__((always_inline)) bblas_dltrsm(const ompssblas_t uplo, const ompssblas_t trans, const ompssblas_t diag, 
		int bm, int bn, int m, int n, double alpha, double *A, int lda, double *B, int ldb) {

	if ( TEST_LOWER(uplo) ) {
		if ( TEST_NTRANSP(trans) ) {
			/* lower, ntrans */
			task_dtrsm(OMPSSBLAS_LEFT, uplo, trans, diag, bm, bn, alpha, &A[0], lda, &B[0], ldb, 1); 

			int i;
			for ( i=bm; i<m; i+=bm ) {
				int dimA = m-i > bm ? bm : m-i;
				task_dgemm(trans, OMPSSBLAS_NTRANSP, dimA, bn, bm, FP_MONE, &A[i], lda, &B[0], ldb, alpha, &B[i], ldb, 1);

				int offs = i + bm * lda; 
				int j;
				for ( j=bm; j<i; j+=bm ) {
					int dimB = m-j > bm ? bm : m-j;
					task_dgemm(trans, OMPSSBLAS_NTRANSP, dimA, bn, dimB, FP_MONE, &A[offs], lda, &B[j], ldb, FP_ONE, &B[i], ldb, 1);
					offs += bm * lda;
				}

				task_dtrsm(OMPSSBLAS_LEFT, uplo, trans, diag, dimA, bn, FP_ONE, &A[offs], lda, &B[i], ldb, 1); 
			}
		} else {
			/* lower, trans */
			int bR = m%bm > 0 ? m%bm : bm;
			int boffs = m - bR;
			int aoffs = boffs * lda + boffs;
			int step = bm;

			task_dtrsm(OMPSSBLAS_LEFT, uplo, trans, diag, bR, bn, alpha, &A[aoffs], lda, &B[boffs], ldb, 1); 

			int j;
			for ( j = boffs-bm; j >= 0; j-=bm, step+=bm ) {
				aoffs -= bm * lda;
				task_dgemm(trans, OMPSSBLAS_NTRANSP, bm, bn, bR, FP_MONE, &A[aoffs], lda, &B[boffs], ldb, alpha, &B[j], ldb, 1);

				int i;
				for ( i = bm; i < step; i+=bm ) {
					task_dgemm(trans, OMPSSBLAS_NTRANSP, bm, bn, bm, FP_MONE, &A[aoffs-i], lda, &B[boffs-i], ldb, FP_ONE, &B[j], ldb, 1);
				}

				task_dtrsm(OMPSSBLAS_LEFT, uplo, trans, diag, bm, bn, FP_ONE, &A[j*lda+j], lda, &B[j], ldb, 1); 
			}
		}
	} else {
		if ( TEST_NTRANSP(trans) ) {
			/* upper, ntrans */
			int bR = m%bm > 0 ? m%bm : bm;
			int boffs = m - bR;
			int aoffs = boffs * lda + boffs;
			int step = bm;

			task_dtrsm(OMPSSBLAS_LEFT, uplo, trans, diag, bR, bn, alpha, &A[aoffs], lda, &B[boffs], ldb, 1);

			int i;
			for ( i = boffs-bm; i >= 0; i-=bm, step+=bm ) {
				aoffs -= bm;
				task_dgemm(trans, OMPSSBLAS_NTRANSP, bm, bn, bR, FP_MONE, &A[aoffs], lda, &B[boffs], ldb, alpha, &B[i], ldb, 1);

				int offs = aoffs;
				int j;
				for ( j = bm; j < step; j+=bm ) {
					offs -= bm * lda;
					task_dgemm(trans, OMPSSBLAS_NTRANSP, bm, bn, bm, FP_MONE, &A[offs], lda, &B[boffs-j], ldb, FP_ONE, &B[i], ldb, 1);
				}
				task_dtrsm(OMPSSBLAS_LEFT, uplo, trans, diag, bm, bn, FP_ONE, &A[i*lda+i], lda, &B[i], ldb, 1);
			}
		} else {
			/* upper, trans */
			task_dtrsm(OMPSSBLAS_LEFT, uplo, trans, diag, bm, bn, alpha, &A[0], lda, &B[0], ldb, 1); 

			int j;
			for ( j=bm; j<m; j+=bm ) {
				int dimA = m-j > bm ? bm : m-j;
				task_dgemm(trans, OMPSSBLAS_NTRANSP, dimA, bn, bm, FP_MONE, &A[j*lda], lda, &B[0], ldb, alpha, &B[j], ldb, 1);

				int i;
				for ( i=bm; i<j; i+=bm ) {
					int dimB = m-i > bm ? bm : m-i;
					task_dgemm(trans, OMPSSBLAS_NTRANSP, dimA, bn, dimB, FP_MONE, &A[j*lda+i], lda, &B[i], ldb, FP_ONE, &B[j], ldb, 1);
				}
				task_dtrsm(OMPSSBLAS_LEFT, uplo, trans, diag, dimA, bn, FP_ONE, &A[j*lda+i], lda, &B[j], ldb, 1); 
			}
		}
	}
}

#endif // __BBLAS_TRSM_H__
