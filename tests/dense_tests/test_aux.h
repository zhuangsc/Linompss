#ifndef __TEST_AUX__
#define __TEST_AUX__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX(a,b) 	((a)>(b)?(a):(b))
#define MIN(a,b) 	((a)<(b)?(a):(b))

#ifdef SINGLE_PRECISION

#define FP			float
#define THRES		1E-5
#define SQRT		sqrtf
#define ABS			fabsf

#define CBLAS_GEMM		cblas_sgemm
#define CBLAS_SYRK		cblas_ssyrk
#define CBLAS_TRSM  	cblas_strsm
#define CBLAS_TRMM		cblas_strmm
#define LAPACK_LANGE	LAPACKE_slange_work
#define LAPACK_GETRF	LAPACKE_sgetrf_work
#define LAPACK_LACPY	LAPACKE_slacpy_work
#define LAPACK_LAMCH	LAPACKE_slamch_work

#define OMPSS_GEMM		ompss_sgemm
#define OMPSS_SYRK		ompss_ssyrk
#define OMPSS_TRSM		ompss_strsm
#define OMPSS_LU		ompss_sgetrf
#define OMPSS_LULL		ompss_sgetrf_ll
#define OMPSS_CHOL		ompss_schol_rl

#endif

#ifdef DOUBLE_PRECISION

#define FP			double
#define THRES		1E-15
#define SQRT		sqrt
#define ABS			fabs

#define CBLAS_GEMM		cblas_dgemm
#define CBLAS_SYRK		cblas_dsyrk
#define CBLAS_TRSM		cblas_dtrsm
#define CBLAS_TRMM		cblas_dtrmm
#define LAPACK_LANGE	LAPACKE_dlange_work
#define LAPACK_GETRF	LAPACKE_dgetrf_work
#define LAPACK_LACPY	LAPACKE_dlacpy_work
#define LAPACK_LAMCH	LAPACKE_dlamch_work

#define OMPSS_GEMM		ompss_dgemm
#define OMPSS_SYRK  	ompss_dsyrk
#define OMPSS_TRSM		ompss_dtrsm
#define OMPSS_LU		ompss_dgetrf
#define OMPSS_LULL		ompss_dgetrf_ll
#define OMPSS_CHOL		ompss_dchol_rl

#endif

static inline void __print_matrix(FP *mat, int m, int n, int lda)
{
	int i;
	for ( i = 0; i < m; ++i ) {
		int j;
		for ( j = 0; j < n; ++j ) {
			printf("%e ", mat[j*lda+i]);
		}
		printf("\n");
	}
	printf("\n\n");
}

static inline FP __vector_2norm(FP *v, int length)
{
	FP x0 = 0;
	int i;
	for( i = 0; i < length; ++i ) {
		x0 += v[i] * v[i];
	}
	return SQRT(x0);
}

static inline FP __attribute__((always_inline)) __mat_1norm(FP *v, int m, int n)
{
	FP maximum = 0.0;
	int j;
	for ( j = 0; j < n; ++j ) {
		FP acc = 0.0;
		int i;
		for ( i = 0; i < m; ++ i) {
			acc += ABS(v[j*m+i]);
		}
		maximum = acc>maximum ? acc : maximum;
	}
	return maximum;

}

static inline FP __attribute__((always_inline)) mat_relerr(FP *v, FP *v0, int m, int n)
{
	int length=m*n;
	FP *vtmp = malloc(length * sizeof(FP));
	int i;
	for ( i = 0; i < length; ++i )
		vtmp[i] = v[i] - v0[i];

	FP Nv0 = __mat_1norm(v0, m, n);
	FP Ntmp = __mat_1norm(vtmp, m, n);
	free(vtmp);
	return Ntmp/Nv0;
}

static inline void __attribute__((always_inline)) GENMAT_SYM_FULL(int m, int n, double *A) 
{
	srand48(time(NULL));

	int j;
	for (j = 0; j < m; ++j ) {
		int i;
		for( i = j; i < n; ++i ) {
			FP dran = drand48();
			A[j*m+i] = A[i*m+j] = dran;
		}
  	}
	for(j = 0; j < m; ++j)
		A[j*m+j] += 10 * m;
}


#endif
