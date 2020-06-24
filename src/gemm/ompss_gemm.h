#ifndef __OMPSS_GEMM_H__
#define __OMPSS_GEMM_H__


#ifdef __cplusplus
extern "C" {
#endif


#ifdef LIBOMPSS_BUILDING
#define LIBOMPSS_DLL_EXPORTED __attribute__((__visibility__("default")))
#else
#define LIBOMPSS_DLL_EXPORTED
#endif


extern LIBOMPSS_DLL_EXPORTED int ompss_hyper_sgemm(int b, int d, int c, int m, int n, int k, float alpha, float *A, int ldimA, float *B, int ldimB, float beta, float *C, int ldimC);

extern LIBOMPSS_DLL_EXPORTED int ompss_hyper_dgemm(int b, int d, int c, int m, int n, int k, float alpha, float *A, int ldimA, float *B, int ldimB, float beta, float *C, int ldimC);

extern LIBOMPSS_DLL_EXPORTED int ompss_sgemm(int b, int d, int c, int m, int n, int k, float alpha, float *A, int ldimA, float *B, int ldimB, float beta, float *C, int ldimC);

extern LIBOMPSS_DLL_EXPORTED int ompss_gemm(int b, int d, int c, int m, int n, int k, float alpha, float *A, int ldimA, float *B, int ldimB, float beta, float *C, int ldimC);


#ifdef __cplusplus
}
#endif


#endif // __OMPSS_GEMM_H__
