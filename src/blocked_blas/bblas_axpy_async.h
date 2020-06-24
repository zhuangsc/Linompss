#ifndef __BBLAS_AXPY_ASYNC_H__
#define __BBLAS_AXPY_ASYNC_H__


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "fptype.h"
#include "task_axpy_async.h"


#ifdef SINGLE_PRECISION

#define BBLAS_AXPY 							bblas_saxpy
#define BBLAS_SCAL_AXPY 					bblas_scal_saxpy
#define BBLAS_SCAL_CPAXPY 					bblas_scal_scpaxpy
#define BBLAS_EXT_AXPY 						bblas_ext_saxpy
#define BBLAS_CPAXPY						bblas_scpaxpy
#define BBLAS_SCAL_CPAXPY_COMB_ASYNC 		bblas_scal_scpaxpy_comb_async
#define BBLAS_SCAL_CPAXPY_COMB_ASYNC_RELEASE 		bblas_scal_scpaxpy_comb_async_release
#define BBLAS_SCAL_CPAXPY_COMB_ASYNC_CONCURRENT		bblas_scal_scpaxpy_comb_async_concurrent
#define BBLAS_CPAXPY_COMB4_ASYNC 			bblas_scpaxpy_comb4_async

#else // DOUBLE_PRECISION

#define BBLAS_AXPY 							bblas_daxpy
#define BBLAS_SCAL_AXPY 					bblas_scal_daxpy
#define BBLAS_SCAL_CPAXPY 					bblas_scal_dcpaxpy
#define BBLAS_EXT_AXPY 						bblas_ext_daxpy
#define BBLAS_CPAXPY						bblas_dcpaxpy
#define BBLAS_SCAL_CPAXPY_COMB_ASYNC 		bblas_scal_dcpaxpy_comb_async
#define BBLAS_SCAL_CPAXPY_COMB_ASYNC_RELEASE 		bblas_scal_dcpaxpy_comb_async_release
#define BBLAS_SCAL_CPAXPY_COMB_ASYNC_CONCURRENT		bblas_scal_dcpaxpy_comb_async_concurrent
#define BBLAS_CPAXPY_COMB4_ASYNC 			bblas_dcpaxpy_comb4_async

#endif


#if 0
/* D is for dummy */
static inline void __attribute__((always_inline)) bblas_saxpy_async(async_t *sync, int bm, int bn, int m, int n, float *Anum, float *Aden, float *X, float *D, float *Y) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int j;
		for ( j=0; j<n; j+=bn ) {
			task_saxpy_async(dotidasync, sync, bm, bn, m, n, &Anum[j], &Aden[j], &X[j*m+i], &D[j*m+i], &Y[j*m+i]);
		}
	}
}
#endif


/* Z = Y + Anum / Aden * X, D is for dummy */
static inline void __attribute__((always_inline)) bblas_scpaxpy_async(int id, async_t *sync, int bm, int bn, int m, int n, float *Anum, float *Aden, float *X, float *D, float *Y, float *Z) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			int idx = j*m+i;
			task_scpaxpy_async(id, sync, c, d, m, n, &Anum[j], &Aden[j], &X[idx], &D[idx], &Y[idx], &Z[idx]);
		}
	}
}


#if 0
static inline void __attribute__((always_inline)) bblas_scal_saxpy_async(async_t *sync, int bm, int bn, int m, int n, float alpha, float *Anum, float *Aden, float *X, float *Y) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int j;
		for ( j=0; j<n; j+=bn ) {
			task_scal_saxpy_async(dotidasync, sync, bm, bn, m, n, alpha, &Anum[j], &Aden[j], &X[j*m+i], &Y[j*m+i]);
		}
	}
}
#endif


static inline void __attribute__((always_inline)) bblas_scal_scpaxpy_async(int id, async_t *sync, int bm, int bn, int m, int n, float alpha, float *Anum, float *Aden, float *X, float *Y, float *Z) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			task_scal_scpaxpy_async(id, sync, c, d, m, n, alpha, &Anum[j], &Aden[j], &X[j*m+i], &Y[j*m+i], &Z[j*m+i]);
		}
	}
}


static inline void __attribute__((always_inline)) bblas_scal_scpaxpy_comb_async(int id, async_t *sync, int p, int bm, int bn, int m, int n, 
	float alpha, float *Anum, float *Aden, float *X1, float *X2, float *Y1, float *Y2, float *Z1, float *Z2) 
{
	int l;
	int i;
	for ( l=0, i=0; i<m; i+=bm, ++l ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			task_scal_scpaxpy_comb_async(id, l, sync, p, c, d, m, n, alpha, &Anum[j], &Aden[j], &X1[j*m+i], &X2[j*m+i], &Y1[j*m+i], &Y2[j*m+i], &Z1[j*m+i], &Z2[j*m+i]);
			//task_scal_scpaxpy_async(dotidasync, sync, c, d, m, n, alpha, &Anum[j], &Aden[j], &X1[j*m+i], &Y1[j*m+i], &Z1[j*m+i]);
		}
	}
}


static inline void __attribute__((always_inline)) bblas_scal_dcpaxpy_comb_async(int id, async_t *sync, int p, int bm, int bn, int m, int n, 
	double alpha, double *Anum, double *Aden, double *X1, double *X2, double *Y1, double *Y2, double *Z1, double *Z2) 
{
	int l;
	int i;
	for ( l=0, i=0; i<m; i+=bm, ++l ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			task_scal_dcpaxpy_comb_async(id, l, sync, p, c, d, m, n, alpha, &Anum[j], &Aden[j], &X1[j*m+i], &X2[j*m+i], &Y1[j*m+i], &Y2[j*m+i], &Z1[j*m+i], &Z2[j*m+i]);
		}
	}
}


static inline void __attribute__((always_inline)) bblas_scal_scpaxpy_comb_async_release(int id, async_t *sync, int p, int bm, int bn, int m, int n, 
	float alpha, float *Anum, float *Aden, float *X1, float *X2, float *Y1, float *Y2, float *Z1, float *Z2) 
{
	int l;
	int i;
	for ( l=0, i=0; i<m; i+=bm, ++l ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			task_scal_scpaxpy_comb_async_release(id, l, sync, p, c, d, m, n, alpha, &Anum[j], &Aden[j], &X1[j*m+i], &X2[j*m+i], &Y1[j*m+i], &Y2[j*m+i], &Z1[j*m+i], &Z2[j*m+i]);
			//task_scal_scpaxpy_async(dotidasync, sync, c, d, m, n, alpha, &Anum[j], &Aden[j], &X1[j*m+i], &Y1[j*m+i], &Z1[j*m+i]);
		}
	}
}


static inline void __attribute__((always_inline)) bblas_scal_dcpaxpy_comb_async_release(int id, async_t *sync, int p, int bm, int bn, int m, int n, 
	double alpha, double *Anum, double *Aden, double *X1, double *X2, double *Y1, double *Y2, double *Z1, double *Z2) 
{
	int l;
	int i;
	for ( l=0, i=0; i<m; i+=bm, ++l ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			task_scal_dcpaxpy_comb_async_release(id, l, sync, p, c, d, m, n, alpha, &Anum[j], &Aden[j], &X1[j*m+i], &X2[j*m+i], &Y1[j*m+i], &Y2[j*m+i], &Z1[j*m+i], &Z2[j*m+i]);
		}
	}
}


static inline void __attribute__((always_inline)) bblas_scal_scpaxpy_comb_async_concurrent(int id, async_t *sync, int p, int bm, int bn, int m, int n, 
	float alpha, float *Anum, float *Aden, float *X1, float *X2, float *Y1, float *Y2, float *Z1, float *Z2) 
{
	int l;
	int i;
	for ( l=0, i=0; i<m; i+=bm, ++l ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			task_scal_scpaxpy_comb_async_concurrent(id, l, sync, p, c, d, m, n, alpha, &Anum[j], &Aden[j], &X1[j*m+i], &X2[j*m+i], &Y1[j*m+i], &Y2[j*m+i], &Z1[j*m+i], &Z2[j*m+i]);
			//task_scal_scpaxpy_async(dotidasync, sync, c, d, m, n, alpha, &Anum[j], &Aden[j], &X1[j*m+i], &Y1[j*m+i], &Z1[j*m+i]);
		}
	}
}


static inline void __attribute__((always_inline)) bblas_scal_dcpaxpy_comb_async_concurrent(int id, async_t *sync, int p, int bm, int bn, int m, int n, 
	double alpha, double *Anum, double *Aden, double *X1, double *X2, double *Y1, double *Y2, double *Z1, double *Z2) 
{
	int l;
	int i;
	for ( l=0, i=0; i<m; i+=bm, ++l ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			task_scal_dcpaxpy_comb_async_concurrent(id, l, sync, p, c, d, m, n, alpha, &Anum[j], &Aden[j], &X1[j*m+i], &X2[j*m+i], &Y1[j*m+i], &Y2[j*m+i], &Z1[j*m+i], &Z2[j*m+i]);
		}
	}
}




#if 0
/* Y = Anum / Aden * X, D is for dummy */
static inline void __attribute__((always_inline)) bblas_daxpy_async(async_t *sync, int bm, int bn, int m, int n, double *Anum, double *Aden, double *X, double *D, double *Y) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int j;
		for ( j=0; j<n; j+=bn ) {
			task_daxpy_async(dotidasync, sync, bm, bn, m, n, &Anum[j], &Aden[j], &X[j*m+i], &D[j*m+i], &Y[j*m+i]);
		}
	}
}
#endif

/* Z = Y + Anum / Aden * X, D is for dummy */
static inline void __attribute__((always_inline)) bblas_dcpaxpy_async(int id, async_t *sync, int bm, int bn, int m, int n, double *Anum, double *Aden, double *X, double *D, double *Y, double *Z) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			int idx = j*m+i;
			task_dcpaxpy_async(id, sync, c, d, m, n, &Anum[j], &Aden[j], &X[idx], &D[idx], &Y[idx], &Z[idx]);
		}
	}
}

#if 0
static inline void __attribute__((always_inline)) bblas_scal_daxpy_async(async_t *sync, int bm, int bn, int m, int n, double alpha, double *Anum, double *Aden, double *X, double *Y) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int j;
		for ( j=0; j<n; j+=bn ) {
			task_scal_daxpy_async(dotidasync, sync, bm, bn, m, n, alpha, &Anum[j], &Aden[j], &X[j*m+i], &Y[j*m+i]);
		}
	}
}

#endif




static inline void __attribute__((always_inline)) bblas_scal_dcpaxpy_async(int id, async_t *sync, int bm, int bn, int m, int n, double alpha, double *Anum, double *Aden, double *X, double *Y, double *Z) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			task_scal_dcpaxpy_async(id, sync, c, d, m, n, alpha, &Anum[j], &Aden[j], &X[j*m+i], &Y[j*m+i], &Z[j*m+i]);
		}
	}
}


static inline void __attribute__((always_inline)) bblas_scpaxpy_comb4_async(int id, async_t *sync, int p, int bm, int bn, int m, int n, float *gamma1, float *gamma2, float *delta, float *sigma2,\
	float *P1, float *V1, float *X1, float *R1, float *S, float *sigma1, float *P2, float *V2, float *X2, float *R2) {

	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int i;
		for ( i=0; i<m; i+=bm ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			int offs = j * m + i;

			task_scpaxpy_comb4_async(id, sync, p, c, d, m, n, gamma1, gamma2, delta, sigma2, &P1[offs], &V1[offs], &X1[offs], &R1[offs], &S[offs], 
				sigma1, &P2[offs], &V2[offs], &X2[offs], &R2[offs]);
		}

		gamma1 += bn;
		gamma2 += bn;
		delta += bn;
		sigma1 += bn;
		sigma2 += bn;
	}
}


static inline void __attribute__((always_inline)) bblas_dcpaxpy_comb4_async(int id, async_t *sync, int p, int bm, int bn, int m, int n, double *gamma1, double *gamma2, double *delta, double *sigma2,\
	double *P1, double *V1, double *X1, double *R1, double *S, double *sigma1, double *P2, double *V2, double *X2, double *R2) {

	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int i;
		for ( i=0; i<m; i+=bm ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			int offs = j * m + i;

			task_dcpaxpy_comb4_async(id, sync, p, c, d, m, n, gamma1, gamma2, delta, sigma2, &P1[offs], &V1[offs], &X1[offs], &R1[offs], &S[offs], 
				sigma1, &P2[offs], &V2[offs], &X2[offs], &R2[offs]);
		}

		gamma1 += bn;
		gamma2 += bn;
		delta += bn;
		sigma1 += bn;
		sigma2 += bn;
	}
}





#endif // __BBLAS_AXPY_ASYNC_H__
