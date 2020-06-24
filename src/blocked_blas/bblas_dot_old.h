#ifndef __BBLAS_DOT_H__
#define __BBLAS_DOT_H__


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fptype.h"
#include "tasks.h"
//#include "check.h"


extern int __ddotids;
extern int __sdotids;
extern int __dnormids;
extern int __snormids;

static inline void __attribute__((always_inline)) dot_sync_print(syncdot_t *tmp){

	printf("Create: %d\n", tmp->create);
	printf("pcnt: %d\n", tmp->pcnt);
	printf("pcompl: %d\n", tmp->pcompl);
	printf("ready: %d\n", tmp->ready);
}

static inline syncdot_t * __attribute__((always_inline)) dot_sync_init(int c, int taskc, void *work) {
	syncdot_t *tmp = (syncdot_t *) work;

	int i;
	for ( i=0; i<c; ++i ) {
		tmp[i].create = -1;
		tmp[i].pcnt = 0;
		tmp[i].pcompl = taskc;
		tmp[i].ready = -1;
		pthread_mutex_init(&tmp[i].mutex, NULL);
	}

	return tmp;
}

static inline void __attribute__((always_inline)) dot_sync_fini(int c, void *work) {
	syncdot_t *tmp = (syncdot_t *) work;

	int i;
	for ( i=0; i<c; ++i ) {
		pthread_mutex_destroy(&tmp[i].mutex);
	}
}


static inline void __attribute__((always_inline)) bblas_sdot(syncdot_t *sync, int bm, int bn, int m, int n, float *X, float *Y, float *result) {
	++__sdotids;
	int j;
	for ( j=0; j<n; j+=bn ) {
		int i;
		for ( i=0; i<m; ) {
			int k = m-i<bm ? m-i : bm;
			task_sdot(__sdotids, sync, k, bn, m, n, &X[j*m+i], &Y[j*m+i], result);
			i += k;
		}
		result += bn;
	}
}


static inline void __attribute__((always_inline)) bblas_ddot(syncdot_t *sync, int bm, int bn, int m, int n, double *X, double *Y, double *result) {
	++__ddotids;
	int j;
	for ( j=0; j<n; j+=bn ) {
		int i;
		for ( i=0; i<m; ) {
			int k = m-i<bm ? m-i : bm;
			task_ddot(__ddotids, sync, k, bn, m, n, &X[j*m+i], &Y[j*m+i], result);
			i += k;
		}
		result += bn;
	}
}


static inline void __attribute__((always_inline)) bblas_s2norm(syncdot_t *sync, int bm, int bn, int m, int n, float *X, float *Y, float *result) {
	++__snormids;
	int j;
	for ( j=0; j<n; j+=bn ) {
		int i;
		for ( i=0; i<m; i+=bm ) {
			task_s2norm(__snormids, sync, bm, bn, m, n, &X[j*m+i], &Y[j*m+i], result);
		}
		result += bn;
	}
}


static inline void __attribute__((always_inline)) bblas_d2norm(syncdot_t *sync, int bm, int bn, int m, int n, double *X, double *Y, double *result) {
	++__dnormids;
	int j;
	for ( j=0; j<n; j+=bn ) {
		int i;
		for ( i=0; i<m; i+=bm ) {
			task_d2norm(__dnormids, sync, bm, bn, m, n, &X[j*m+i], &Y[j*m+i], result);
		}
		result += bn;
	}
}


#endif // __BBLAS_DOT_H__
