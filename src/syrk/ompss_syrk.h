#ifndef __OMPSSSYRK_H__
#define __OMPSSSYRK_H__

#ifdef __cplusplus
extern "C" {
#endif


#include "fpblas.h"
#include "fpmatr.h"


//#ifdef LIBOMPSS_BUILDING
#define LIBOMPSS_DLL_EXPORTED __attribute__((__visibility__("default")))
//#else
//#define LIBOMPSS_DLL_EXPORTED
//#endif


extern LIBOMPSS_DLL_EXPORTED void ompss_ssyrk(ompssblas_t uplo, ompssblas_t trans, int b, int n, int k, float alpha, float *A, int lda, float beta, float *C, int ldc);

extern LIBOMPSS_DLL_EXPORTED void ompss_dsyrk(ompssblas_t uplo, ompssblas_t trans, int b, int n, int k, double alpha, double *A, int lda, double beta, double *C, int ldc);


#ifdef __cplusplus
}
#endif

#endif //__OMPSSSYRK_H__
