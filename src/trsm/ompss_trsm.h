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


extern LIBOMPSS_DLL_EXPORTED void ompss_strsm(ompssblas_t side, ompssblas_t uplo, ompssblas_t trans, ompssblas_t diag, \
		int m, int n, int b, float alpha, float *A, int lda, float *B, int ldb);

extern LIBOMPSS_DLL_EXPORTED void ompss_dtrsm(ompssblas_t side, ompssblas_t uplo, ompssblas_t trans, ompssblas_t diag, \
		int m, int n, int b, double alpha, double *A, int lda, double *B, int ldb);


#ifdef __cplusplus
}
#endif

#endif //__OMPSSSYRK_H__
