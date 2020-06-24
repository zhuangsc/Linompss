#ifndef __CSRMMB_H__
#define __CSRMMB_H__


#include "hb.h"


#ifdef DOUBLE_PRECISION

#define CSRMMB	dcsrmmb

#else

#define CSRMMB	scsrmmb

#endif


int scsrmmb(int b, int c, int d, int n, float alpha, hbmat_t *Ahbh, float *B, int ldb, float beta, float *C, int ldc);
int dcsrmmb(int b, int c, int d, int n, double alpha, hbmat_t *Ahbh, double *B, int ldb, double beta, double *C, int ldc);


#endif // __CSRMMB_H__
