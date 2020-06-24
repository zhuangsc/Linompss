#ifndef __GENMAT_H__
#define __GENMAT_H__


#ifdef SINGLE_PRECISION

#define GENMAT					sgenmat
#define GENMAT_ZERO				sgenmat_zero
#define GENMAT_IP				sgenmat_ip
#define GENMAT_SPD				sgenmat_spd
#define GENMAT_HYPER_SPD		sgenmat_hyper_spd
#define GENMAT_SYM_FULL			sgenmat_sym_full
#define GENMAT_SYM_LTRIANG		sgenmat_sym_ltriang
#define GENMAT_EYE				sgenmat_eye
#define GENMAT_LTRIANG_ONES		sgenmat_ltriang_ones
#define GENMAT_SYMDD			sgenmat_symdd
#define GENMAT_COVARIANT		sgenmat_covariant
#define GENMAT_COSTAS			sgenmat_costas

#else

#define GENMAT					dgenmat
#define GENMAT_ZERO				dgenmat_zero
#define GENMAT_IP				dgenmat_ip
#define GENMAT_SPD				dgenmat_spd
#define GENMAT_HYPER_SPD		dgenmat_hyper_spd
#define GENMAT_SYM_FULL			dgenmat_sym_full
#define GENMAT_SYM_LTRIANG		dgenmat_sym_ltriang
#define GENMAT_EYE				dgenmat_eye
#define GENMAT_LTRIANG_ONES		dgenmat_ltriang_ones
#define GENMAT_SYMDD			dgenmat_symdd
#define GENMAT_COVARIANT		dgenmat_covariant
#define GENMAT_COSTAS			dgenmat_costas

#endif


#define LIBBBLAS_EXPORT __attribute__((__visibility__("default")))


#include <stdlib.h>


extern LIBBBLAS_EXPORT float * sgenmat_spd(int n, int ldim); 
extern LIBBBLAS_EXPORT double * dgenmat_spd(int n, int ldim); 

extern LIBBBLAS_EXPORT void sgenmat_hyper_spd(int n, int tn, int nleft, int b, float *A); 
extern LIBBBLAS_EXPORT void dgenmat_hyper_spd(int n, int tn, int nleft, int b, double *A); 

extern LIBBBLAS_EXPORT float* sgenmat(int m, int n, int ld);
extern LIBBBLAS_EXPORT double* dgenmat(int m, int n, int ld);

extern LIBBBLAS_EXPORT void sgenmat_ip(float *A, int m, int n, int scale);
extern LIBBBLAS_EXPORT void dgenmat_ip(double *A, int m, int n, int scale);

extern LIBBBLAS_EXPORT void sgenmat_sym_full(int m, float *A);
extern LIBBBLAS_EXPORT void dgenmat_sym_full(int m, double *A);

extern LIBBBLAS_EXPORT void sgenmat_sym_ltriang(int m, float *A);
extern LIBBBLAS_EXPORT void dgenmat_sym_ltriang(int m, double *A);

extern LIBBBLAS_EXPORT void* dgenmat_eye(int m, int ldim);
extern LIBBBLAS_EXPORT void* sgenmat_eye(int m, int ldim);

extern LIBBBLAS_EXPORT void* dgenmat_zero(int m, int ldim);
extern LIBBBLAS_EXPORT void* sgenmat_zero(int m, int ldim);

extern LIBBBLAS_EXPORT void sgenmat_ltriang_ones(int m, float *A);
extern LIBBBLAS_EXPORT void dgenmat_ltriang_ones(int m, double *A);

extern LIBBBLAS_EXPORT float* sgenmat_covariant(size_t n, float p, float k);
extern LIBBBLAS_EXPORT double* dgenmat_covariant(size_t n, double p, double k);

extern LIBBBLAS_EXPORT float* sgenmat_costas(size_t m, size_t n);
extern LIBBBLAS_EXPORT double* dgenmat_costas(size_t m, size_t n);


#endif // __GENMAT_H__
