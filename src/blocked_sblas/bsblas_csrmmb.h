#ifndef __BSBLAS_CSRMMB_H__
#define __BSBLAS_CSRMMB_H__


#include "fptype.h"
#include "hb.h"
#include "task_csrmmb.h"
#include "matfprint.h"
#include "async_struct.h"


#ifdef SINGLE_PRECISION

#define BSBLAS_CSRMMB				bsblas_scsrmmb
#define BSBLAS_CSRMMB_PRED			bsblas_scsrmmb_pred

#else

#define BSBLAS_CSRMMB				bsblas_dcsrmmb
#define BSBLAS_CSRMMB_PRED			bsblas_dcsrmmb_pred

#endif


/* A hbh csr, B and C dense col major */
static inline void __attribute__((always_inline)) bsblas_scsrmmb(int b, int n, float alpha, hbmat_t *Ahbh, float *B, int ldb, float beta, float *C, int ldc) 
{
	hbmat_t *A = Ahbh->orig;
	int m = A->m;
	int k = A->n;

	int M = Ahbh->m;
	int N = (n + b - 1 ) / b;
	int K = Ahbh->n;

	int *vptr = Ahbh->vptr;
	int *vpos = Ahbh->vpos;
	hbmat_t **vval = Ahbh->vval;
	int offs = vptr[0] == 0 ? 0 : 1;

	/* change F to C for row-major, adjust leading dimensions */
	int cmaj = 1;

	int I;
	for ( I = 0; I < M; ++I ) {
		int Le = vptr[I+1] - offs;
		int L0 = vptr[I] - offs;

		int j;
		for ( j = 0; j < n; j+=b ) {
			int left = n - j;
			int jb = b > left? left : b;

			int k = vpos[L0] - offs;
			int boffs = cmaj ? j * ldb + k * b : k * b * ldb + j;
			int coffs = cmaj ? j * ldc + I * b : I * b * ldc + j;

			float *bC = &(C[coffs]);
			float *bB = &(B[boffs]);

			task_scsrmmb(jb, alpha, "GLNF", vval[L0], bB, ldb, beta, bC, ldc, 0, NULL);

			int L;
			for ( L = L0+1; L < Le; ++L ) {
				int k = vpos[L] - offs;
				int boffs = cmaj ? j * ldb + k * b : k * b * ldb + j;
				int coffs = cmaj ? j * ldc + I * b : I * b * ldc + j;

				float *bB = &(B[boffs]);

				task_scsrmmb(jb, alpha, "GLNF", vval[L], bB, ldb, FP_ONE, bC, ldc, 0, NULL);
			}
		}
	}
}


static inline void __attribute__((always_inline)) bsblas_dcsrmmb(int b, int n, double alpha, hbmat_t *Ahbh, double *B, int ldb, double beta, double *C, int ldc) 
{
	hbmat_t *A = Ahbh->orig;
	int m = A->m;
	int k = A->n;

	int M = Ahbh->m;
	int N = (n + b - 1 ) / b;
	int K = Ahbh->n;

	int *vptr = Ahbh->vptr;
	int *vpos = Ahbh->vpos;
	hbmat_t **vval = Ahbh->vval;
	int offs = vptr[0] == 0 ? 0 : 1;

	/* change F to C for row-major, adjust leading dimensions */
	int cmaj = 1;

	int I;
	for ( I = 0; I < M; ++I ) {
		int Le = vptr[I+1] - offs;
		int L0 = vptr[I] - offs;

		int j;
		for ( j = 0; j < n; j+=b ) {
			int left = n - j;
			int jb = b > left? left : b;

			int k = vpos[L0] - offs;
			int boffs = cmaj ? j * ldb + k * b : k * b * ldb + j;
			int coffs = cmaj ? j * ldc + I * b : I * b * ldc + j;

			double *bC = &(C[coffs]);
			double *bB = &(B[boffs]);

			task_dcsrmmb(jb, alpha, "GLNF", vval[L0], bB, ldb, beta, bC, ldc, 0, NULL);

			int L;
			for ( L = L0+1; L < Le; ++L ) {
				int k = vpos[L] - offs;
				int boffs = cmaj ? j * ldb + k * b : k * b * ldb + j;
				int coffs = cmaj ? j * ldc + I * b : I * b * ldc + j;

				double *bB = &(B[boffs]);

				task_dcsrmmb(jb, alpha, "GLNF", vval[L], bB, ldb, FP_ONE, bC, ldc, 0, NULL);
			}
		}
	}
}


static inline async_stat_t __attribute__((always_inline)) bsblas_scsrmmb_pred(int b, int n, float alpha, hbmat_t *Ahbh, float *B, int ldb, float beta, float *C, int ldc, \
	int it, int est, int prev, volatile int *cond, double **alpha1, double tol) 
{
 	async_stat_t status = STAT_AHEAD;

	hbmat_t *A = Ahbh->orig;
	int m = A->m;
	int k = A->n;

	int M = Ahbh->m;
	int N = (n + b - 1 ) / b;
	int K = Ahbh->n;

	int *vptr = Ahbh->vptr;
	int *vpos = Ahbh->vpos;
	hbmat_t **vval = Ahbh->vval;
	int offs = vptr[0] == 0 ? 0 : 1;

	/* change F to C for row-major, adjust leading dimensions */
	int cmaj = 1;

	int I;
	for ( I = 0; I < M; ++I ) {
		int Le = vptr[I+1] - offs;
		int L0 = vptr[I] - offs;

		int j;
		for ( j = 0; j < n; j+=b ) {
			int left = n - j;
			int jb = b > left? left : b;

			int k = vpos[L0] - offs;
			int boffs = cmaj ? j * ldb + k * b : k * b * ldb + j;
			int coffs = cmaj ? j * ldc + I * b : I * b * ldc + j;

			float *bC = &(C[coffs]);
			float *bB = &(B[boffs]);

			if ( status == STAT_AHEAD ) {
//                status = ASMAN_BREAK(it, est, prev, alpha1, cond, 0, tol);
				printf("csrmmb error!\n");
                if ( status == STAT_BROKEN || status == STAT_CONVERGED ) {
                    return status;
                }
            }


			task_scsrmmb(jb, alpha, "GLNF", vval[L0], bB, ldb, beta, bC, ldc, 0, NULL);

			int L;
			for ( L = L0+1; L < Le; ++L ) {
				int k = vpos[L] - offs;
				int boffs = cmaj ? j * ldb + k * b : k * b * ldb + j;
				int coffs = cmaj ? j * ldc + I * b : I * b * ldc + j;

				float *bB = &(B[boffs]);

				task_scsrmmb(jb, alpha, "GLNF", vval[L], bB, ldb, FP_ONE, bC, ldc, 0, NULL);
			}
		}
	}

	return status;
}


static inline async_stat_t __attribute__((always_inline)) bsblas_dcsrmmb_pred(int b, int n, double alpha, hbmat_t *Ahbh, double *B, int ldb, double beta, double *C, int ldc, \
	int it, int est, int prev, volatile int *cond, double **alpha1, double tol) 
{
 	async_stat_t status = STAT_AHEAD;

	hbmat_t *A = Ahbh->orig;
	int m = A->m;
	int k = A->n;

	int M = Ahbh->m;
	int N = (n + b - 1 ) / b;
	int K = Ahbh->n;

	int *vptr = Ahbh->vptr;
	int *vpos = Ahbh->vpos;
	hbmat_t **vval = Ahbh->vval;
	int offs = vptr[0] == 0 ? 0 : 1;

	/* change F to C for row-major, adjust leading dimensions */
	int cmaj = 1;

	int I;
	for ( I = 0; I < M; ++I ) {
		int Le = vptr[I+1] - offs;
		int L0 = vptr[I] - offs;

		int j;
		for ( j = 0; j < n; j+=b ) {
			int left = n - j;
			int jb = b > left? left : b;

			int k = vpos[L0] - offs;
			int boffs = cmaj ? j * ldb + k * b : k * b * ldb + j;
			int coffs = cmaj ? j * ldc + I * b : I * b * ldc + j;

			double *bC = &(C[coffs]);
			double *bB = &(B[boffs]);

			if ( status == STAT_AHEAD ) {
//                status = ASMAN_BREAK(it, est, prev, alpha1, cond, 0, tol);
				printf("csrmmb error!\n");
                if ( status == STAT_BROKEN || status == STAT_CONVERGED ) {
                    return status;
                }
            }


			task_dcsrmmb(jb, alpha, "GLNF", vval[L0], bB, ldb, beta, bC, ldc, 0, NULL);

			int L;
			for ( L = L0+1; L < Le; ++L ) {
				int k = vpos[L] - offs;
				int boffs = cmaj ? j * ldb + k * b : k * b * ldb + j;
				int coffs = cmaj ? j * ldc + I * b : I * b * ldc + j;

				double *bB = &(B[boffs]);

				task_dcsrmmb(jb, alpha, "GLNF", vval[L], bB, ldb, FP_ONE, bC, ldc, 0, NULL);
			}
		}
	}

	return status;
}



#endif // __BSBLAS_CSRMMB_H__


