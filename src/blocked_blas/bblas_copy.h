#ifndef __BBLAS_COPY_H__
#define __BBLAS_COPY_H__


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fptype.h"
#include "task_copy.h"
#include "selfsched.h"
#include "async_struct.h"


#ifdef SINGLE_PRECISION

#define BBLAS_COPY 				bblas_scopy
#define BBLAS_COPY_SCHED 		bblas_scopy_sched

#else // DOUBLE_PRECISION

#define BBLAS_COPY 				bblas_dcopy
#define BBLAS_COPY_SCHED 		bblas_dcopy_sched

#endif


static inline void __attribute__((always_inline)) bblas_scopy(int p, int bm, int bn, int m, int n, float *X, float *Y) 
{
	int i;
	for ( i=0; i<m; i+=bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			task_scopy(p, c, d, m, n, &X[j*m+i], &Y[j*m+i]);
		}
	}
}


static inline void __attribute__((always_inline)) bblas_dcopy(int p, int bm, int bn, int m, int n, double *X, double *Y) 
{
	int i;
	for ( i=0; i<m; i+=bm ) {
		int cs = m - i;
		int c = cs < bm ? cs : bm;

		int j;
		for ( j=0; j<n; j+=bn ) {
			int ds = n - j;
			int d = ds < bn ? ds : bn;

			task_dcopy(p, c, d, m, n, &X[j*m+i], &Y[j*m+i]);
		}
	}
}

#endif
