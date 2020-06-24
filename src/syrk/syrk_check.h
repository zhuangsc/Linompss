#ifndef __SYRK_CHECK_H__
#define __SYRK_CHECK_H__


#include "fptype.h"
#include "fpblas.h"
#include "fpmatr.h"
#include "densutil.h"

fp_t syrk_check(int check, ompssblas_t uplo, ompssblas_t trans, int b, int n, int k, fp_t alpha, fp_t *A, int lda, fp_t beta, fp_t *C, int ldc, fp_t *Cc);


#endif // __SYRK_CHECK_H__
