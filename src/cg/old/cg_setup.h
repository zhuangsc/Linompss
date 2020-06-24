#ifndef __CG_SETUP_H__
#define __CG_SETUP_H__


#include "fptype.h"


int cg_setup(int n, int bm, int bn, int s, void **A, fp_t **x, fp_t **xstar, fp_t **rhs, fp_t **work);  
int cg_cleanup(int n, int s, void *A, fp_t *x, fp_t *xstar, fp_t *rhs, fp_t *work);


#endif // __CG_SETUP_H__
