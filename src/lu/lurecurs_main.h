#ifndef __LURECURS_MAIN_H__
#define  __LURECURS_MAIN_H__

#include "fptype.h"

extern int IPIVi;
extern fp_t *Astart;

static inline __attribute__((always_inline)) void lu_recurs_start(fp_t *A) {
	IPIVi = 0;
	Astart = A;
}


void lu_recurs(int depth,int skip, int ldim, int n, int t, fp_t *A, int *IPIV);

void lu_recurs_l(int ldim, fp_t *A, int *IPIV);


#endif // __LURECURS_MAIN_H__
