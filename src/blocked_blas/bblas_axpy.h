#ifndef __BBLAS_AXPY_H__
#define __BBLAS_AXPY_H__


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "fptype.h"
#include "task_axpy.h"


#ifdef SINGLE_PRECISION

#define BBLAS_AXPY 					bblas_saxpy
#define BBLAS_AXPY3 				bblas_saxpy3
#define BBLAS_SCAL_AXPY 			bblas_scal_saxpy
#define BBLAS_SCAL_CPAXPY 			bblas_scal_scpaxpy
#define BBLAS_EXT_AXPY 				bblas_ext_saxpy
#define BBLAS_EXTM_AXPY 			bblas_extm_saxpy
#define BBLAS_CPAXPY				bblas_scpaxpy
#define BBLAS_CPAXPY_COMB			bblas_scpaxpy_comb
#define BBLAS_CPAXPY_COMB4          bblas_scpaxpy_comb4

#else // DOUBLE_PRECISION

#define BBLAS_AXPY 					bblas_daxpy
#define BBLAS_AXPY3 				bblas_daxpy3
#define BBLAS_SCAL_AXPY 			bblas_scal_daxpy
#define BBLAS_SCAL_CPAXPY 			bblas_scal_dcpaxpy
#define BBLAS_EXT_AXPY 				bblas_ext_daxpy
#define BBLAS_EXTM_AXPY 			bblas_extm_daxpy
#define BBLAS_CPAXPY				bblas_dcpaxpy
#define BBLAS_CPAXPY_COMB			bblas_dcpaxpy_comb
#define BBLAS_CPAXPY_COMB4          bblas_dcpaxpy_comb4

#endif


#if 0
/* Z = Y - SA * X */
static inline void __attribute__((always_inline)) bblas_ext_daxpy(int bm, int bn, int m, int n, double SA, double *X, double *Y, double *Z) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int j;
		for ( j=0; j<n; j+=bn) {
			task_ext_daxpy(bm, bn, m, n, SA, &X[j*m+i], &Y[j*m+i], &Z[j*m+i]);
		}
	}
}
#endif


static inline void __attribute__((always_inline)) bblas_extm_daxpy(int p, int bm, int bn, int m, int n, double *SAnum, double *SAden, double *X, double *Y, double *Z) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			task_extm_daxpy(c, d, m, n, &SAnum[j], &SAden[j], &X[j*m+i], &Y[j*m+i], &Z[j*m+i], p);
		}
	}
}


/* Z = Y + Anum / Aden * X */
static inline void __attribute__((always_inline)) bblas_scpaxpy(int bm, int bn, int m, int n, float *Anum, float *Aden, float *X, float *Y, float *Z) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			int idx = j*m+i;
			task_scpaxpy(c, d, m, n, &Anum[j], &Aden[j], &X[idx], &Y[idx], &Z[idx]);
		}
	}
}


/* Z = Y + Anum / Aden * X */
static inline void __attribute__((always_inline)) bblas_dcpaxpy(int bm, int bn, int m, int n, double *Anum, double *Aden, double *X, double *Y, double *Z) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			int idx = j*m+i;
			task_dcpaxpy(c, d, m, n, &Anum[j], &Aden[j], &X[idx], &Y[idx], &Z[idx]);
		}
	}
}




#if 0
static inline void __attribute__((always_inline)) bblas_saxpy(int bm, int bn, int m, int n, float *Anum, float *Aden, float *X, float *Y) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int j;
		for ( j=0; j<n; j+=bn ) {
			task_saxpy(bm, bn, m, n, &Anum[j], &Aden[j], &X[j*m+i], &Y[j*m+i]);
		}
	}
}

static inline void __attribute__((always_inline)) bblas_scal_saxpy(int bm, int bn, int m, int n, float alpha, float *Anum, float *Aden, float *X, float *Y) {
        int i;
        for ( i=0; i<m; i+=bm ) {
                int j;
                for ( j=0; j<n; j+=bn ) {
                        task_scal_saxpy(bm, bn, m, n, alpha, &Anum[j], &Aden[j], &X[j*m+i], &Y[j*m+i]);
                }
        }
}
#endif


/* Z = Y + Anum / Aden * X */
static inline void __attribute__((always_inline)) bblas_scal_scpaxpy(int bm, int bn, int m, int n, float alpha, float *Anum, float *Aden, float *X, float *Y, float *Z) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			task_scal_scpaxpy(c, d, m, n, alpha, &Anum[j], &Aden[j], &X[j*m+i], &Y[j*m+i], &Z[j*m+i]);
		}
	}
}


#if 0
static inline void __attribute__((always_inline)) bblas_daxpy(int bm, int bn, int m, int n, double *Anum, double *Aden, double *X, double *Y) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int j;
		for ( j=0; j<n; j+=bn ) {
			task_daxpy(bm, bn, m, n, &Anum[j], &Aden[j], &X[j*m+i], &Y[j*m+i]);
		}
	}
}
#endif


static inline void __attribute__((always_inline)) bblas_sdaxpy(int bm, int bn, int m, int n, float a, float *x, double *y) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			task_sdaxpy(c, d, m, n, a, &x[j*m+i], &y[j*m+i]);
		}
	}
}


static inline void __attribute__((always_inline)) bblas_sdcpaxpy(int bm, int bn, int m, int n, double a, double *x, double *y, double *z) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			task_sdcpaxpy(c, d, m, n, a, &x[j*m+i], &y[j*m+i], &z[j*m+i]);
		}
	}
}


static inline void __attribute__((always_inline)) bblas_scal_dcpaxpy(int bm, int bn, int m, int n, double alpha, double *Anum, double *Aden, double *X, double *Y, double *Z) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			task_scal_dcpaxpy(c, d, m, n, alpha, &Anum[j], &Aden[j], &X[j*m+i], &Y[j*m+i], &Z[j*m+i]);
		}
	}
}


#if 0
/* Z = Y - SA * X */
static inline void __attribute__((always_inline)) bblas_ext_saxpy(int bm, int bn, int m, int n, float SA, float *X, float *Y, float *Z) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int j;
		for ( j=0; j<n; j+=bn) {
			task_ext_saxpy(bm, bn, m, n, SA, &X[j*m+i], &Y[j*m+i], &Z[j*m+i]);
		}
	}
}
#endif


static inline void __attribute__((always_inline)) bblas_extm_saxpy(int p, int bm, int bn, int m, int n, float *SAnum, float *SAden, float *X, float *Y, float *Z) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			task_extm_saxpy(c, d, m, n, &SAnum[j], &SAden[j], &X[j*m+i], &Y[j*m+i], &Z[j*m+i], p);
		}
	}
}

#if 0
static inline void __attribute__((always_inline)) bblas_extm_saxpy_sched(int z, int p, int bm, int bn, int m, int n, float *SAnum, float *SAden, selfsched_t *sched, float *X, float *Y, float *Z) 
{
	int l;
	int i;
	for ( i=0, l=0; i<m; i+=bm, ++l ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			task_extm_saxpy_sched(z, p, c, d, m, n, &SAnum[j], &SAden[j], l, sched, &X[j*m+i], &Y[j*m+i], &Z[j*m+i]);
		}
	}
}
#endif




static inline void __attribute__((always_inline)) bblas_scpaxpy_comb(int bm, int bn, int m, int n, float alpha, float *Anum, float *Aden, float *X1, float *X2, 
	float *Y1, float *Y2, float *Z1, float *Z2) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			task_scpaxpy_comb(c, d, m, n, alpha, &Anum[j], &Aden[j], &X1[j*m+i], &X2[j*m+i], &Y1[j*m+i], &Y2[j*m+i], &Z1[j*m+i], &Z2[j*m+i]);
		}
	}
}

static inline void __attribute__((always_inline)) bblas_dcpaxpy_comb(int bm, int bn, int m, int n, double alpha, double *Anum, double *Aden, double *X1, double *X2, 
	double *Y1, double *Y2, double *Z1, double *Z2) {
	int i;
	for ( i=0; i<m; i+=bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			task_dcpaxpy_comb(c, d, m, n, alpha, &Anum[j], &Aden[j], &X1[j*m+i], &X2[j*m+i], &Y1[j*m+i], &Y2[j*m+i], &Z1[j*m+i], &Z2[j*m+i]);
		}
	}
}


static inline void __attribute__((always_inline)) bblas_scpaxpy_comb4(int bm, int bn, int m, int n, float *gamma1, float *gamma2, float *delta, float *sigma2,\
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

            task_scpaxpy_comb4(c, d, m, n, gamma1, gamma2, delta, sigma2, &P1[offs], &V1[offs], &X1[offs], &R1[offs], &S[offs],
                sigma1, &P2[offs], &V2[offs], &X2[offs], &R2[offs]);
        }

        gamma1 += bn;
        gamma2 += bn;
        delta += bn;
        sigma1 += bn;
        sigma2 += bn;
    }
}


static inline void __attribute__((always_inline)) bblas_saxpy3(int bm, int bn, int m, int n, float *gamman, float *gammad, float *alpha, 
	float *z, float *s, float *p, float *q, float *w, float *xp, float *x, float *rp, float *r) 
{
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int i;
		for ( i=0; i<m; i+=bm ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			int offs = j * m + i;

			task_saxpy4(c, d, m, n, gamman, gammad, alpha, &z[offs], &s[offs], &p[offs], &q[offs], &w[offs], &xp[offs], &x[offs], &rp[offs], &r[offs]);
		}

		gamman += bn;
		gammad += bn;
		alpha += bn;
	}
}

static inline void __attribute__((always_inline)) bblas_daxpy3(int bm, int bn, int m, int n, double *gamman, double *gammad, double *alpha, 
	double *z, double *s, double *p, double *q, double *w, double *xp, double *x, double *rp, double *r) 
{
	int j;
	for ( j=0; j<n; j+=bn ) {
		int ds = n - j;
		int d = ds < bn ? ds : bn;

		int i;
		for ( i=0; i<m; i+=bm ) {
			int cs = m - i;
			int c = cs < bm ? cs : bm;

			int offs = j * m + i;

			task_daxpy4(c, d, m, n, gamman, gammad, alpha, &z[offs], &s[offs], &p[offs], &q[offs], &w[offs], &xp[offs], &x[offs], &rp[offs], &r[offs]);
		}

		gamman += bn;
		gammad += bn;
		alpha += bn;
	}
}



static inline void __attribute__((always_inline)) bblas_dcpaxpy_comb4(int bm, int bn, int m, int n, double *gamma1, double *gamma2, double *delta, double *sigma2,\
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

            task_dcpaxpy_comb4(c, d, m, n, gamma1, gamma2, delta, sigma2, &P1[offs], &V1[offs], &X1[offs], &R1[offs], &S[offs],
                sigma1, &P2[offs], &V2[offs], &X2[offs], &R2[offs]);
        }

        gamma1 += bn;
        gamma2 += bn;
        delta += bn;
        sigma1 += bn;
        sigma2 += bn;
    }
}


#endif // __BBLAS_AXPY_H__
