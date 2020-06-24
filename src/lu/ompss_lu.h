#ifndef __OMPSSLU_H__
#define __OMPSSLU_H__

#ifdef __cplusplus
extern "C" {
#endif


#include "fpblas.h"
//#include "fpmatr.h"


#define LIBOMPSS_DLL_EXPORTED __attribute__((__visibility__("default")))

extern LIBOMPSS_DLL_EXPORTED void ompss_sgetrf(int m, int n, int b, float *A, int lda, int *IPIV);

extern LIBOMPSS_DLL_EXPORTED void ompss_dgetrf(int m, int n, int b, double *A, int lda, int *IPIV);

extern LIBOMPSS_DLL_EXPORTED void ompss_sgetrf_ll(int m, int n, int b, float *A, int lda, int *IPIV);

extern LIBOMPSS_DLL_EXPORTED void ompss_dgetrf_ll(int m, int n, int b, double *A, int lda, int *IPIV);

//extern LIBOMPSS_DLL_EXPORTED void ompss_sgetrf2(int m, int n, int b, float *A, int lda, int *IPIV);
//extern LIBOMPSS_DLL_EXPORTED void ompss_dgetrf2(int m, int n, int b, double *A, int lda, int *IPIV);


#ifdef __cplusplus
}
#endif

#endif //__OMPSSSYRK_H__
