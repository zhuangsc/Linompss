#ifndef __CSRMM_H__
#define __CSRMM_H__


#include "hb.h"
#include "selfsched.h"


#ifdef DOUBLE_PRECISION

#define CSRMM	dcsrmm

#else

#define CSRMM	scsrmm

#endif


int scsrmm(int b, int c, int d, int n, float alpha, hbmat_t *Ahbh, float *B, int ldb, float beta, float *C, int ldc, selfsched_t *sched);
int dcsrmm(int b, int c, int d, int n, double alpha, hbmat_t *Ahbh, double *B, int ldb, double beta, double *C, int ldc, selfsched_t *sched);


#endif // __CSRMM_H__
