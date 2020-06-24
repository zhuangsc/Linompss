#ifndef __ITREF_SETUP_H__
#define __ITREF_SETUP_H__


#include "fptype.h"


int itref_setup(int n, int bm, int bn, int s, void **A, double *x[], double **rhs, float **work);  
int itref_cleanup(int n, int s, void *A, double *x[], double *rhs, float *work);


#endif // __ITREF_SETUP_H__
