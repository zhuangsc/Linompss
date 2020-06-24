#ifndef __MM_H__
#define __MM_H__


#include "fptype.h"


#ifdef DOUBLE_PRECISION
#define MM	dmm
#else
#define MM	smm
#endif


int dmm(int p, int b, int c, int d, int m, int k, int n, int transa, int transb, double *A, double *B, double *C);
int smm(int p, int b, int c, int d, int m, int k, int n, int transa, int transb, float *A, float *B, float *C);


#endif // __MM_H__
