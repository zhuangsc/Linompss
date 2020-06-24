#ifndef __AS_MAN__
#define __AS_MAN__


#include "async_struct.h"


#include <stdlib.h>
#include <stdio.h>


typedef struct {
	int istar;
	int itop;
	int iprev;
	int *istack;
} asman_t;
/* Keep track of 4 not necessarily different solutions: 
	1) the best approximation 
	2) the previous approx 
	3) the one under computation
	4) the next one */

static inline asman_t* __attribute__((always_inline)) asman_init(int c) {
	asman_t *a = malloc(sizeof(asman_t));

	int *st = a->istack = malloc(c * sizeof(int));
	int i;
	for ( i=0; i<c; ++i) {
		st[i] = i;
	}

	a->istar = -1;
	a->itop = 0;
	a->iprev = -1;
	
	return a;
}

static inline void __attribute__((always_inline)) asman_fini(asman_t *a) {
	free( a->istack );
	free(a);
	a = NULL;
}


static inline void __attribute__((always_inline)) asman_update(asman_t *a, int star, int prev) {
    a->istar = star;
    a->iprev = prev;
}

static inline int __attribute__((always_inline)) asman_nexti(asman_t *a) {
	return a->istack[a->itop++];
}


static inline int __attribute__((always_inline)) asman_best(asman_t *a) {
	return a->istar;
}


#if 0
static inline void __attribute__((always_inline)) asman_repair() {
	printf("repair %i %i %i %i\n", istack[0], istack[1], istack[2], istack[3]);
	printf("star %i top %i prev %i\n", istar, itop, iprev);
}
#endif

#if SINGLE_PRECISION
#define ASMAN_BREAK 	asman_sbreak
#else
#define ASMAN_BREAK 	asman_dbreak
#endif


extern async_stat_t asman_sbreak(int it, int est/*, int prev */, asman_t *asman, void **alpha1, async_t *sync, int force, double tol, float *residuals);
extern async_stat_t asman_dbreak(int it, int est/*, int prev */, asman_t *asman, void **alpha1, async_t *sync, int force, double tol, double *residuals);


#endif // __AS_MAN__
