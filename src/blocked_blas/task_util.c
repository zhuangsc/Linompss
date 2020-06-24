#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "fptype.h"
#include "task_util.h"


#ifdef DOUBLE_PRECISION

#define __t_clear		task_dclear
#define __t_set			task_dset

#else

#define __t_clear		task_sclear
#define __t_set			task_sset

#endif


#define DRYRUN	0



void __t_clear(int n, fp_t *a) {
#if ! DRYRUN
	int i;
	for ( i=0; i<n; ++i) {
		a[i] = 0.0;
	}
#endif
}


void __t_set(int bm, int bn, int m, int n, fp_t v, fp_t *A) {
#if ! DRYRUN
	int i;
	for ( i=0; i<bm; ++i) {
		int j;
		for ( j=0; j<bn; ++j ) {
			A[j*m+i] = v;
		}
	}
#endif
}

