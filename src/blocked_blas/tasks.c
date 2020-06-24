#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fptype.h"
#include "fpblas.h"
#include "tasks.h"
#include "blas.h"



#ifdef DOUBLE_PRECISION

#define __t_copy 		task_dcopy
#define __t_gemm			task_dgemm
#define __t_gemm0beta	task_dgemm0beta
#define __t_clear		task_dclear
#define __t_div			task_ddiv
#define __t_ext_div		task_ext_ddiv
#define __t_scal			task_dscal
#define __t_inv_scal			task_inv_dscal
#define __t_nrm2			task_dnrm2
#define __t_set			task_dset

int dlltest;

#else

#define __t_copy 		task_scopy
#define __t_gemm			task_sgemm
#define __t_gemm0beta	task_sgemm0beta
#define __t_clear		task_sclear
#define __t_div			task_sdiv
#define __t_ext_div		task_ext_sdiv
#define __t_scal			task_sscal
#define __t_inv_scal			task_inv_sscal
#define __t_nrm2			task_snrm2
#define __t_set			task_sset

#endif

#define DRYRUN	0


void __t_copy(int bm, int bn, int m, int n, fp_t *X, fp_t *Y) {
#if ! DRYRUN
	//verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, X);
	int j;
	for ( j=0; j<bn; ++j ) {
		BLAS_cp(bm, &X[j*m], I_ONE, &Y[j*m], I_ONE);
	}
	//verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, Y);
#endif
}


void __t_gemm0beta(int bm, int bn, int bk, fp_t alpha, fp_t *A, int m, fp_t *B, int n, fp_t *C, int k) {
#if ! DRYRUN
	//verify_2df_array(__FILE__, __LINE__, M, K, BSM, BSK, A);
	//verify_2df_array(__FILE__, __LINE__, K, N, BSK, BSN, B);
	//verify_2df_array(__FILE__, __LINE__, M, N, BSM, BSN, C);
	BLAS_gemm(MAT_NOTRANSP, MAT_NOTRANSP, bm, bn, bk, alpha, A, m, B, n, FP_NOUGHT, C, k);
	//verify_2df_array(__FILE__, __LINE__, M, N, BSM, BSN, C);
#endif
}




void __t_clear(int n, fp_t *a) {
#if ! DRYRUN
	int i;
	for ( i=0; i<n; ++i) {
		a[i] = 0.0;
	}
#endif
}



void __t_div(int bn, fp_t *X, fp_t *Y) {
#if ! DRYRUN
	int i;
	for ( i=0; i<bn; ++i) {
		Y[i] /= X[i];
	}
	
	#if CORRECT_DIVS
		_Bool warned = 0;
		for (int i=0; i < BS; i++) {
			if (isnan(Y[i]) || isinf(Y[i])) {
				Y[i] = 0.0;
				if (!warned) {
					#ifdef CELL_SUPERSCALAR_H
						fprintf(stdout, "Fixed div value in task %i\n", getCurrentTaskId());
					#else
						fprintf(stdout, "Fixed div value\n");
					#endif
					warned = 1;
				}
			}
		}
	#endif
#endif
}


void __t_ext_div(int bn, fp_t *X, fp_t *Y, fp_t *Z) {
#if ! DRYRUN
	int i;
	for ( i=0; i<bn; ++i) {
		Z[i] = X[i] / Y[i];
	}
#endif
}


void __t_scal(int bm, int bn, int m, int n, fp_t *SA, fp_t *X) {
#if ! DRYRUN
	//verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, X);
	//verify_2df_array(__FILE__, __LINE__, 1, N, 1, BSB, SA);
	int j;
	for ( j=0; j<bn; ++j ) {
		printf("scaling with %f\n", SA[j]);
		BLAS_scal(bm, SA[j], X, I_ONE);
		X += m;
	}
	//verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, X);
#endif
}

void __t_inv_scal(int bm, int bn, int m, int n, fp_t *SA, fp_t *X) {
#if ! DRYRUN
	//verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, X);
	//verify_2df_array(__FILE__, __LINE__, 1, N, 1, BSB, SA);
	int j;
	for ( j=0; j<bn; ++j ) {
		fp_t fact = 1 / SA[j];
		BLAS_scal(bm, fact, X, I_ONE);
		X += m;
	}
	//verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, X);
#endif
}



void __t_set(int bm, int bn, int m, int n, fp_t v, fp_t *A) {
#if ! DRYRUN
	int i;
	for ( i=0; i<bm; ++i) {
		int j;
		for ( j=0; j<bn; ++j ) {
			A[j*m+i] = v;
		}
	}
#endif
}

void __t_nrm2(int bm, int bn, fp_t *A, int ldim, fp_t *norm2) {
#if 0
	int j;
	for ( j=0; j<bn; ++j ) {
    		norm2[j] = BLAS_nrm2(&bm, A, &one);
		A += m;
	}
#endif
}


#if 0

void dcopy_2dtr_tile(int BSA, int BSB, int M, int N, double *X, double *Y) {
	verify_2dd_array(__FILE__, __LINE__, M, N, BSA, BSB, X);
	for (int i=0; i<BSA; i++) {
		dcopy_(&BSB, &X[i*N], &one, &Y[i*N], &one);
	}
	verify_2dd_array(__FILE__, __LINE__, M, N, BSA, BSB, Y);
}



void dscal_2dtr_tile(int BSA, int BSB, int M, int N, double DA, double *X) {
	verify_2dd_array(__FILE__, __LINE__, M, N, BSA, BSB, X);
	for (int i=0; i<BSA; i++) {
		dscal_(&BSB, &DA, &X[i*N], &one);
	}
	verify_2dd_array(__FILE__, __LINE__, M, N, BSA, BSB, X);
}


void sscal_2dtr_tile(int BSA, int BSB, int M, int N, float A, float *X) {
	verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, X);
	for (int i=0; i<BSA; i++) {
		sscal_(&BSB, &A, &X[i*N], &one);
	}
	verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, X);
}


void saxpy_2dtrsa_tile(int BSA, int BSB, int M, int N, float SA, float *X, float *Y) {
	verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, X);
	verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, Y);
	for (int i=0; i<BSA; i++) {
		saxpy_(&BSB, &SA, &X[i*N], &one, &Y[i*N], &one);
	}
	verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, Y);
}


void saxpy_2dtr_ext_tile(int BSA, int BSB, int M, int N, float *SA, float *X, float *Y, float *Z) {
	verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, X);
	verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, Y);
	verify_2df_array(__FILE__, __LINE__, 1, N, 1, BSB, SA);
	for (int j=0; j<BSB; j++) {
		scopy_(&BSA, &Y[j], &N, &Z[j], &N);
		saxpy_(&BSA, &SA[j], &X[j], &N, &Z[j], &N);
	}
	verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, Z);
}




void dset_2dtr_tile(int BSA, int BSB, int M, int N, double value, double *a) {
	for (int i=0; i<BSA; i++) {
		for (int j=0; j<BSB; j++) {
			a[j+i*N] = value;
		}
	}
	verify_2dd_array(__FILE__, __LINE__, M, N, BSA, BSB, a);
}


void sset_2dtr_tile(int BSA, int BSB, int M, int N, float value, float *a) {
	for (int i=0; i<BSA; i++) {
		for (int j=0; j<BSB; j++) {
			a[j+i*N] = value;
		}
	}
	verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, a);
}


void d2s_2dtr_tile(int BSA, int BSB, int M, int N, double *X, float *Y) {
	verify_2dd_array(__FILE__, __LINE__, M, N, BSA, BSB, X);
	for (int i=0; i < BSA; i++) {
		for (int j=0; j < BSB; j++) {
			Y[j+i*N] = (float)X[j+i*N];
		}
	}
	verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, Y);
}


void s2d_2dtr_tile(int BSA, int BSB, int M, int N, float *X, double *Y) {
	verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, X);
	for (int i=0; i < BSA; i++) {
		for (int j=0; j < BSB; j++) {
			Y[j+i*N] = (double)X[j+i*N];
		}
	}
	verify_2dd_array(__FILE__, __LINE__, M, N, BSA, BSB, Y);
}


void dps_2dtr_tile(int BSA, int BSB, int M, int N, float *X, double *Y) {
	verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, X);
	verify_2dd_array(__FILE__, __LINE__, M, N, BSA, BSB, Y);
	for (int i=0; i < BSA; i++) {
		for (int j=0; j < BSB; j++) {
			Y[j+i*N] += (double)X[j+i*N];
		}
	}
	verify_2dd_array(__FILE__, __LINE__, M, N, BSA, BSB, Y);
}


void init_tile(int BSA, int BSB, int M, int N, int i, int j, double *A) {
	for (int I=0; I<BSA; I++) {
		for (int J=0; J<BSB; J++) {
			if (I+i != J+j) {
				double tmpval = (double)(I+i-J-j);
				A[J+I*N] = 1.0/(tmpval*tmpval*tmpval*tmpval);
			} else {
				A[J+I*N] = 1.0 + sqrt((double)(I+i+1));
			}
		}
	}
	verify_2dd_array(__FILE__, __LINE__, M, N, BSA, BSB, A);
}


void dnrm2_tr_task(int N, int s, double *R, double *result) {
	verify_2dd_array(__FILE__, __LINE__, N, s, N, 1, R);
	*result = dnrm2_(&N, R, &s);
}
#endif




