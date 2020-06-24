#include "chol_kernels.h"

#include <stdio.h>

#include "fptype.h"
#include "fpblas.h"
#include "fplapack.h"
#include "blas.h"


#if 1
void GEMM_TASK( int b, int t, fp_t *A, int lda, fp_t *B, int ldb, fp_t *C, int ldc) {
#if 1
	BLAS_gemm(MAT_NOTRANSP, MAT_TRANSP, b, b, b, FP_MONE, A, lda, B, ldb, FP_ONE, C, ldc);
#else
	int tm = t*ldm;
	int k;
	for ( k=0; k<b; k+=t ) {
		int kb = min(b-k,t);

		GEMM("No transpose", "Transpose", &b, &kb, &b, &mone, A, &ldm, B, &ldm, &one, C, &ldm);

		C += tm;
		B += t;
	}
#endif
}
#endif
			      	

void SYRK_TASK(int b, fp_t *A, int lda, fp_t *C, int ldc) {
	BLAS_syrk(TRIANG_LOWER, MAT_NOTRANSP, b, b, FP_MONE, A, lda, FP_ONE, C, ldc);
}


void POTRF_TASK(int b, int t, fp_t *A, int ldm) {
#if 1
	int info;
	LAPACK_potrf(LAPACK_TRIANG_LOWER, b, A, ldm, info);
	if ( info!=0 ) {
		fprintf(stderr,"error in DPOTRF (%i)\n",info);
	}
#else
	fp_t dmone = -1.00;
	fp_t done = 1.00;
	int tm = t * ldm;

	int r = b; 

	int k;
	for ( k =0; k < b; k+=t ) {
		int kb = min(b-k, t);
		fp_t *B = A + kb;
		r -= kb;

		//printf("dpotrf b %i\n", kb);
		int info;
		POTRF("Lower", &kb, A, &ldm, &info);
		if (info!=0) {
			fprintf(stderr,"error: DPOTRF (%i)\n", info);
		}

		r = r>0 ? r : 0;
		//printf("dtrsm %i %i\n", r, kb);
		TRSM("Right", "Lower", "Transpose", "Non-unit", &r, &kb, &done, A, &ldm, B, &ldm);
		//print_matrix_blocked( stdout, "A", m, m, m, m, m, t, Aorig );

		fp_t *Bj = B;
		fp_t *Cj = B + tm; 
		int jr = r;
		int j;
		for ( j = k+t; j < b; j+=t ) {
			int jb = min(b-j, t);
			jr -= jb;
			jr = jr > 0? jr : 0;

			fp_t *Dj = Cj + jb;
			fp_t *Ej = Bj + jb;

			//printf("dsyrk %i %i\n",jb,kb);
			SYRK("Lower", "No Transpose", &jb, &kb, &dmone, Bj, &ldm, &done, Cj, &ldm);
			Cj += tm + t; 
	
			//printf("dgemm %i %i %i\n", jr, jb, jb);
			GEMM("No transpose", "Transpose", &jr, &jb, &t, &dmone, Ej, &ldm, Bj, &ldm, &done, Dj, &ldm);
			//print_matrix_blocked( stdout, "A", m, m, m, m, m, t, Aorig );

			Bj += t;
		}

		A = B + tm;
	}
#endif
}
			      	

void TRSM_TASK( int b, int t, fp_t *A, fp_t *B, int ldm) {
#if 1
	BLAS_trsm(SIDE_RIGHT, TRIANG_LOWER, MAT_TRANSP, DIAG_NUNIT, b, b, FP_ONE, A, ldm, B, ldm);
#else
	fp_t dmone = -1.00;
	int tm = t*ldm;
	int km = b;
	int k; 
	for ( k=0; k<b; k+=t ) {
		int kb = min(b-k,t);
		
		TRSM("Right", "Lower", "Transpose", "Non-unit", &b, &kb, &one, A, &ldm, B, &ldm);

		fp_t *C = B + tm;
		km -= t;
		A += t;

		GEMM("No Transpose", "Transpose", &b, &km, &t, &dmone, B, &ldm, A, &ldm, &one, C, &ldm);

		B = C ;
		A += tm;
	}
#endif
}
