#ifndef __BBLAS_DOT_H__
#define __BBLAS_DOT_H__


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fptype.h"
#include "task_dot.h"


#ifdef SINGLE_PRECISION

#define BBLAS_DOT 			bblas_sdot
#define BBLAS_DOT_PURE		bblas_sdot_pure
#define BBLAS_DOT3 			bblas_sdot3
#define BBLAS_DOT4 			bblas_sdot4

#else // DOUBLE_PRECISION

#define BBLAS_DOT 			bblas_ddot
#define BBLAS_DOT_PURE		bblas_ddot_pure
#define BBLAS_DOT3 			bblas_ddot3
#define BBLAS_DOT4 			bblas_ddot4

#endif


/* FIX: does not work if n > bn */
static inline void __attribute__((always_inline)) bblas_sdot(int id, async_t *sync, int bm, int bn, int m, int n, float *X, float *Y, float *result) 
{
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int idx;
		int i;
		for ( i=0, idx=0; i<m; i+=bm, ++idx ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			task_sdot(id, idx, sync, c, d, m, n, &X[j*m+i], &Y[j*m+i], result);
		}
		result += bn;
	}
}


static inline void __attribute__((always_inline)) bblas_ddot(int id, async_t *sync, int bm, int bn, int m, int n, double *X, double *Y, double *result) 
{
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int idx;
		int i;
		for ( i=0, idx=0; i<m; i+=bm, ++idx ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			task_ddot(id, idx, sync, c, d, m, n, &X[j*m+i], &Y[j*m+i], result);
		}
		result += bn;
	}
}


static inline void __attribute__((always_inline)) bblas_sdot_pure(int p, int bm, int bn, int m, int n, float *X, float *Y, float *result) 
{
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int idx;
		int i;
		for ( i=0, idx=0; i<m; i+=bm, ++idx ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			task_sdot_pure(p, c, d, m, n, &X[j*m+i], &Y[j*m+i], result);
		}
		result += bn;
	}
}


static inline void __attribute__((always_inline)) bblas_ddot_pure(int p, int bm, int bn, int m, int n, double *X, double *Y, double *result) 
{
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int idx;
		int i;
		for ( i=0, idx=0; i<m; i+=bm, ++idx ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			task_ddot_pure(p, c, d, m, n, &X[j*m+i], &Y[j*m+i], result);
		}
		result += bn;
	}
}


/* r1 = < X, X > and r2 = < X, Y > */
static inline void __attribute__((always_inline)) bblas_sdot3(int id, async_t *sync, int bm, int bn, int m, int n, float *X, float *Y, float *r1, float *r2) {
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int idx;
		int i;
		for ( i=0, idx=0; i<m; i+=bm, ++idx ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			task_sdot3(id, idx, sync, c, d, m, n, &X[j*m+i], &Y[j*m+i], r1, r2);
		}
		r1 += bn;
		r2 += bn;
	}
}

/* r1 = < X, X > and r2 = < X, Y > */
static inline void __attribute__((always_inline)) bblas_ddot3(int id, async_t *sync, int bm, int bn, int m, int n, double *X, double *Y, double *r1, double *r2) {
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int idx;
		int i;
		for ( i=0, idx=0; i<m; i+=bm, ++idx ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			task_ddot3(id, idx, sync, c, d, m, n, &X[j*m+i], &Y[j*m+i], r1, r2);
		}
		r1 += bn;
		r2 += bn;
	}
}


/* r1 = < X, X > and r2 = < X, Y > */
static inline void __attribute__((always_inline)) bblas_sdot4(int id, async_t *sync, int bm, int bn, int m, int n, float *prevr1, float *prevalpha, 
	float *X, float *Y, float *r1, float *r2, float *alpha) {
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int i;
		for ( i=0; i<m; i+=bm ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			task_sdot4(id, sync, c, d, m, n, prevr1, prevalpha, &X[j*m+i], &Y[j*m+i], r1, r2, alpha);
		}
		r1 += bn;
		r2 += bn;
	}
}

static inline void __attribute__((always_inline)) bblas_ddot4(int id, async_t *sync, int bm, int bn, int m, int n, double *prevr1, double *prevalpha, 
	double *X, double *Y, double *r1, double *r2, double *alpha) {
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int i;
		for ( i=0; i<m; i+=bm ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			task_ddot4(id, sync, c, d, m, n, prevr1, prevalpha, &X[j*m+i], &Y[j*m+i], r1, r2, alpha);
		}
		r1 += bn;
		r2 += bn;
	}
}


#endif
