#ifndef __ARRAY_H__
#define __ARRAY_H__

#include "../fp/fptype.h"

static inline void __attribute__((always_inline)) array_d2s(hbmat_t* A, fp_t* peel, int col) {
	int* vptr = A->vptr;
	int* vpos = A->vpos;
	fp_t* vval = A->vval;
	int i;
	for ( i = vptr[col]; i < vptr[col+1]; i++ ) {
		vval[i] = peel[vpos[i]];
	}
}


static inline void __attribute__((always_inline)) array_s2d(hbmat_t* A, fp_t* peel, int col) {
	int* vptr = A->vptr;
	int* vpos = A->vpos;
	fp_t* vval = A->vval;
	int i;
	for ( i = vptr[col]; i < vptr[col+1]; i++ ) {
		peel[vpos[i]] = vval[i];
	}
}


static inline void __attribute__((always_inline)) array_clear(fp_t* peel, int elemc) {
	int i;
	for ( i=0; i<elemc; i++ ) {
		peel[i] = 0;
	}
}


static inline void __attribute__((always_inline)) array_clear2(hbmat_t *A, fp_t* peel, int col) {
	int *vptr = A->vptr;	
	int *vpos = A->vpos;	

	int i;
	for ( i=vptr[col]; i<vptr[col+1]; ++i ) {
		peel[vpos[i]] = 0;
	}
}




#endif // __ARRAY_H__
