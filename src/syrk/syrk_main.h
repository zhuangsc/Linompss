#ifndef __SYRK_MAIN_H__
#define __SYRK_MAIN_H__

#include "fptype.h"
#include "fpmatr.h"
#include "fpblas.h"


#ifdef SINGLE_PRECISION

#define SYRK_MAIN		ssyrk_main

#endif 

#ifdef DOUBLE_PRECISION

#define SYRK_MAIN		dsyrk_main

#endif


void dsyrk_main(ompssblas_t uplo, ompssblas_t trans, int b, int n, int k, double alpha, double *A, int lda, double beta, double *C, int ldc);
void ssyrk_main(ompssblas_t uplo, ompssblas_t trans, int b, int n, int k, float  alpha, float  *A, int lda, float  beta, float  *C, int ldc);


#endif // __SYRK_MAIN_H__
