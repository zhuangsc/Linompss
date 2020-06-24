#ifndef __BLOCKED_BLAS_H__
#define __BLOCKED_BLAS_H__


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fptype.h"
#include "task_util.h"


static inline void __attribute__((always_inline)) bblas_sclear(int bn, int n, float *a) 
{
	int i;
	for ( i=0; i<n; i+=bn ) {
		task_sclear(bn, &a[i]);
	}
}


static inline void __attribute__((always_inline)) bblas_sset(int bm, int bn, int m, int n, float v, float *A) 
{
	int i;
	for ( i=0; i<m; i+=bm ) {
		int mr = m-i < bm ? m-i : bm;
		int j;
		for ( j=0; j<n; j+=bn ) {
			int nr = n-j < bn ? n-j : bn;
			task_sset(mr, nr, m, n, v, &A[j*m+i]);
		}
	}
}


static inline void __attribute__((always_inline)) bblas_dclear(int bn, int n, double *a) 
{
	int i;
	for ( i=0; i<n; i+=bn ) {
		task_dclear(bn, &a[i]);
	}
}


static inline void __attribute__((always_inline)) bblas_dset(int bm, int bn, int m, int n, double v, double *A) 
{
	int i;
	for ( i=0; i<m; i+=bm ) {
		int mr = m-i < bm ? m-i : bm;
		int j;
		for ( j=0; j<n; j+=bn ) {
			int nr = n-j < bn ? n-j : bn;
			task_dset(mr, nr, m, n, v, &A[j*m+i]);
		}
	}
}


#endif
