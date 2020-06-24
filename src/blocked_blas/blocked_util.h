#ifndef __BLOCKED_UTIL_H__
#define __BLOCKED_UTIL_H__


#include "task_convert.h"


static inline void __attribute__((always_inline)) block_d2s(int bm, int bn, int m, int n, double *D, float *S) {
	int i;
	for ( i=0; i<m; i+=bm) {
		int j;
		for ( j=0; j<n; j+=bn ) {
		    task_d2s(bm, bn, m, n, &D[j*m+i], &S[j*m+i]);
		}
    	}
}


static inline void __attribute__((always_inline)) block_s2d(int bm, int bn, int m, int n,  float *S, double *D) {
	int i;
	for ( i=0; i<m; i+=bm) {
		int j;
		for ( j=0; j<n; j+=bn ) {
    			task_s2d(bm, bn, m, n, &S[j*m+i], &D[j*m+i]);
        	}
    	}
}


#endif // __BLOCKED_UTIL_H__
