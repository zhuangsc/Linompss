#ifndef __BBLAS_GEMV_H__
#define __BBLAS_GEMV_H__


#include "fptype.h"
#include "tasks.h"


/* y = alpha * A * x + beta * y */
/* y: mx1 , A: mxn, x: nx1 */
static inline void __attribute__((always_inline)) bblas_sgemv(int bm, int bn, int m, int n, float alpha, float *A, float *x, float beta, float *y) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int xi = 0;
		task_sgemv(bm, bn, m, n, alpha, &A[0*m+i], &x[xi], beta, &y[i]);
		int j;
		for ( j=bn; j<n; j+=bn ) {
			xi += bn;
			task_sgemv(bm, bn, m, n, alpha, &A[j*m+i], &x[xi], FP_ONE, &y[i]);
		}
	}
}


static inline void __attribute__((always_inline)) bblas_dgemv(int bm, int bn, int m, int n, double alpha, double *A, double *x, double beta, double *y) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int xi = 0;
		task_dgemv(bm, bn, m, n, alpha, &A[0*m+i], &x[xi], beta, &y[i]);
		int j;
		for ( j=bn; j<n; j+=bn ) {
			xi += bn;
			task_dgemv(bm, bn, m, n, alpha, &A[j*m+i], &x[xi], FP_ONE, &y[i]);
		}
	}
}


#endif // __BBLAS_GEMV_H__
