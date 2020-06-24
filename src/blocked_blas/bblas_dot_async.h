#ifndef __BBLAS_ASYNC_H__
#define __BBLAS_ASYNC_H__


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fptype.h"
#include "task_dot_async.h"
#include "async_struct.h"
//#include "log.h"


#ifdef SINGLE_PRECISION

#define BBLAS_DOT_ASYNC 			bblas_sdot_async
#define BBLAS_DOT_SCHED_ASYNC 		bblas_sdot_sched_async
#define BBLAS_DOT_SCHED_ASYNC_RELEASE 		bblas_sdot_sched_async_release
#define BBLAS_DOT_SCHED_ASYNC_CONCURRENT	bblas_sdot_sched_async_concurrent
#define BBLAS_DOT3_ASYNC 			bblas_sdot3_async
#define BBLAS_DOT3_ASYNC 			bblas_sdot3_async

#else // DOUBLE_PRECISION

#define BBLAS_DOT_ASYNC 			bblas_ddot_async
#define BBLAS_DOT_SCHED_ASYNC 		bblas_ddot_sched_async
#define BBLAS_DOT_SCHED_ASYNC_RELEASE 		bblas_ddot_sched_async_release
#define BBLAS_DOT_SCHED_ASYNC_CONCURRENT	bblas_ddot_sched_async_concurrent
#define BBLAS_DOT3_ASYNC 			bblas_ddot3_async

#endif


//extern int dotidasync;


static inline void __attribute__((always_inline)) bblas_sdot_async(int id, async_t *sync, int p, int bm, int bn, int m, int n, float *X, float *Y, float *result) {
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int i;
		for ( i=0; i<m; i+=bm ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			task_sdot_async(id, sync, p, c, d, m, n, &X[j*m+i], &Y[j*m+i], result);
		}
		result += bn;
	}
}

static inline void __attribute__((always_inline)) bblas_ddot_async(int id, async_t *sync, int p, int bm, int bn, int m, int n, double *X, double *Y, double *result) {
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int i;
		for ( i=0; i<m; i+=bm ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			task_ddot_async(id, sync, p, c, d, m, n, &X[j*m+i], &Y[j*m+i], result);
		}
		result += bn;
	}
}


static inline void __attribute__((always_inline)) bblas_sdot_sched_async(int id, async_t *sync, int p, int bm, int bn, int m, int n, float *X, float *Y, float *result, int release) 
{
	//++dotidasync;
	
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int idx;
		int i;
		for ( i=0, idx=0; i<m; i+=bm, ++idx ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			task_sdot_sched_async(id, idx, sync, p, c, d, m, n, &X[j*m+i], &Y[j*m+i], result, release);
		}
		result += bn;
	}
}

static inline void __attribute__((always_inline)) bblas_ddot_sched_async(int id, async_t *sync, int p, int bm, int bn, int m, int n, double *X, double *Y, double *result, int release) 
{
	//++dotidasync;
	
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int idx;
		int i;
		for ( i=0, idx=0; i<m; i+=bm, ++idx ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			task_ddot_sched_async(id, idx, sync, p, c, d, m, n, &X[j*m+i], &Y[j*m+i], result, release);
		}
		result += bn;
	}
}

static inline void __attribute__((always_inline)) bblas_sdot_sched_async_release(int id, async_t *sync, int p, int bm, int bn, int m, int n, float *X, float *Y, float *result, int release) 
{
	//++dotidasync;
	
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int idx;
		int i;
		for ( i=0, idx=0; i<m; i+=bm, ++idx ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			task_sdot_sched_async_release(id, idx, sync, p, c, d, m, n, &X[j*m+i], &Y[j*m+i], result, release);
		}
		result += bn;
	}
}

static inline void __attribute__((always_inline)) bblas_ddot_sched_async_release(int id, async_t *sync, int p, int bm, int bn, int m, int n, double *X, double *Y, double *result, int release) 
{
	//++dotidasync;
	
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int idx;
		int i;
		for ( i=0, idx=0; i<m; i+=bm, ++idx ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			task_ddot_sched_async_release(id, idx, sync, p, c, d, m, n, &X[j*m+i], &Y[j*m+i], result, release);
		}
		result += bn;
	}
}


static inline void __attribute__((always_inline)) bblas_sdot_sched_async_concurrent(int id, async_t *sync, int p, int bm, int bn, int m, int n, float *X, float *Y, float *result, int release, int *bitmap) 
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

			task_sdot_sched_async_concurrent(id, idx, sync, p, c, d, m, n, &X[j*m+i], &Y[j*m+i], result, release, bitmap);
		}
		result += bn;
	}
}

static inline void __attribute__((always_inline)) bblas_ddot_sched_async_concurrent(int id, async_t *sync, int p, int bm, int bn, int m, int n, double *X, double *Y, double *result, int release, int *bitmap) 
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

			task_ddot_sched_async_concurrent(id, idx, sync, p, c, d, m, n, &X[j*m+i], &Y[j*m+i], result, release, bitmap);
		}
		result += bn;
	}
}



/* r1 = < X, X > and r2 = < X, Y > */
static inline void __attribute__((always_inline)) bblas_sdot3_async(int id, async_t *sync, int p, int bm, int bn, int m, int n, float *X, float *Y, float *r1, float *r2, int release) {
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int idx;
		int i;
		for ( i=0, idx=0; i<m; i+=bm, ++idx ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			task_sdot3_async(id, idx, sync, p, c, d, m, n, &X[j*m+i], &Y[j*m+i], r1, r2, release);
		}
		r1 += bn;
		r2 += bn;
	}
}


static inline void __attribute__((always_inline)) bblas_ddot3_async(int id, async_t *sync, int p, int bm, int bn, int m, int n, double *X, double *Y, double *r1, double *r2, int release) {
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int idx;
		int i;
		for ( i=0, idx=0; i<m; i+=bm, ++idx ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			task_ddot3_async(id, idx, sync, p, c, d, m, n, &X[j*m+i], &Y[j*m+i], r1, r2, release);
		}
		r1 += bn;
		r2 += bn;
	}
}


#endif
