#ifndef __OMPSSMM_H__
#define __OMPSSMM_H__

#ifdef __cplusplus
extern "C" {
#endif


#include "hb.h"
#include "fpblas.h"
#include "fpmatr.h"


//#ifdef LIBOMPSS_BUILDING
#define LIBOMPSS_DLL_EXPORTED __attribute__((__visibility__("default")))
//#else
//#define LIBOMPSS_DLL_EXPORTED
//#endif


/* 	A: (in)		SPD matrix in hypermatrix form
	L: (out)	lower-triangular factor */
extern LIBOMPSS_DLL_EXPORTED int ompss_sgemm(ompssblas_t transa, ompssblas_t transb, int b, int d, int c, int m, int n, int k, float alpha, const float *A, int lda, const float *B, int ldb, float beta, float *C, int ldc);

/* 	A: (in)		SPD matrix in hypermatrix form
	L: (out)	lower-triangular factor */
extern LIBOMPSS_DLL_EXPORTED int ompss_dgemm(ompssblas_t transa, ompssblas_t transb, int b, int d, int c, int m, int n, int k, double alpha, const double *A, int lda, const double *B, int ldb, double beta, double *C, int ldc);


extern LIBOMPSS_DLL_EXPORTED int ompss_dcsrmb(int b, int d, int c, int n, double alpha, hbmat_t *Acsr, const double *B, int ldb, double beta, double *C, int ldc);


extern LIBOMPSS_DLL_EXPORTED int ompss_scsrmb(int b, int d, int c, int n, float alpha, hbmat_t *Acsr, const float *B, int ldb, float beta, float *C, int ldc);


#if 0
extern LIBOMPSS_DLL_EXPORTED int ompss_sgemm(int b, int d, int c, int m, int n, int k, float alpha, float *A, int ldimA, float *B, int ldimB, float beta, float *C, int ldimC);

extern LIBOMPSS_DLL_EXPORTED int ompss_gemm(int b, int d, int c, int m, int n, int k, float alpha, float *A, int ldimA, float *B, int ldimB, float beta, float *C, int ldimC);
#endif


#ifdef __cplusplus
}
#endif

#endif //__OMPSSMM_H__
