#ifndef __BLAS_SYRK_H__
#define __BLAS_SYRK_H__


/* C = alpha * A * A**T + beta * C, with C symmetric */
static inline __attribute__((always_inline)) bblas_dsyrk(fpblasd_t uplo, fpblasd_t transa, int b, int n, int k, double alpha, const double *A, int lda, 
double beta, double *C, int ldc) {
	int nmb = n - b;
	int ps = m * b;

	int i;
	for ( i=0; i<n; i+=b ) {
		int jstop = i - b;
		int j;
		for ( j=0; j<jstop; j+=b ) {
			double *Cb = &C[j*m+i];

			double *Ar = &A[i];
			double *Ac = &A[j*ps];

			task_dgemm(MAT_NOTRANSP, MAT_TRANSP, b, b, b, alpha, Ar, lda, Ac, lda, beta, Cb, ldc);

			int l;
			for ( l=b; l<k; l+=b ) {
				Ar += ps;
				Ac += b;

				task_dgemm(MAT_NOTRANSP, MAT_TRANSP, b, b, b, alpha, Ar, lda, Ac, ldb, FP_ONE, Cb, ldc);
			}
		}

		double *Cd = &C[j*m+j];
		double *Ad = &A[i];

		task_dsyrk(uplo, MAT_NOTRANSP, b, b, alpha, Ad, lda, beta, Cd, ldc);

		int l;
		for ( l=b; l<k; l+=b ) {
			Ad += ps;

			task_dsyrk(uplo, MAT_NOTRANSP, b, b, alpha, Ad, lda, FP_ONE, Cd, ldc);
		}
	}
}

	
#endif // __BLAS_SYRK_H__
