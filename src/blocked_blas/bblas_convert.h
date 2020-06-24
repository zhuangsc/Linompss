#ifndef __BBLAS_CONVERT_H__
#define __BBLAS_CONVERT_H__


#include "task_convert.h"

static inline void __attribute__((always_inline)) bblas_d2s(int bm, int bn, int m, int n, double *D, float *S) 
{
	int i;
	for ( i=0; i<m; i+=bm) {
		int mr = m-i < bm ? m-i : bm;
		int j;
		for ( j=0; j<n; j+=bn ) {
			int nr = n-j < bn ? n-j : bn;
		    task_d2s(mr, nr, m, n, &D[j*m+i], &S[j*m+i]);
		}
	}
}


static inline void __attribute__((always_inline)) bblas_s2d(int bm, int bn, int m, int n,  float *S, double *D) 
{
	int i;
	for ( i=0; i<m; i+=bm) {
		int mr = m-i < bm ? m-i : bm;
		int j;
		for ( j=0; j<n; j+=bn ) {
			int nr = n-j < bn ? n-j : bn;
			task_s2d(mr, nr, m, n, &S[j*m+i], &D[j*m+i]);
		}
	}
}

#endif // __BBLAS_CONVERT_H__
