#ifndef __MATFPRINT_H__
#define __MATFPRINT_H__


#include "hb.h"


#ifdef SINGLE_PRECISION

#define FPRINT_CSR2MM			fprint_scsr2mm
#define FPRINT_CSC2MM			fprint_scsc2mm
#define FPRINT_DENSE2MM			fprint_sdense2mm
#define PRINT_DENSE2MM			print_sdense2mm
#define FPRINT_SUFF_DENSE2MM	fprint_suff_sdense2mm
#define PRINT_HB				print_shb

#else

#define FPRINT_CSR2MM			fprint_dcsr2mm
#define FPRINT_CSC2MM			fprint_dcsc2mm
#define FPRINT_DENSE2MM			fprint_ddense2mm
#define PRINT_DENSE2MM			print_ddense2mm
#define FPRINT_SUFF_DENSE2MM	fprint_suff_ddense2mm
#define PRINT_HB				print_dhb

#endif


#define LIBSBBLAS_EXPORT __attribute__((__visibility__("default")))


extern LIBSBBLAS_EXPORT void fprint_scsr2mm(const char *fname, int m, const int *vptr, const int *vpos, const float *vval); 
extern LIBSBBLAS_EXPORT void fprint_dcsr2mm(const char *fname, int m, const int *vptr, const int *vpos, const double *vval);

extern LIBSBBLAS_EXPORT void fprint_scsc2mm(const char *fname, int m, const int *vptr, const int *vpos, const float *vval);
extern LIBSBBLAS_EXPORT void fprint_dcsc2mm(const char *fname, int m, const int *vptr, const int *vpos, const double *vval);

extern LIBSBBLAS_EXPORT void fprint_sdense2mm(const char *fname, const char *name, int m, int n, const float *A, int lda);
extern LIBSBBLAS_EXPORT void fprint_ddense2mm(const char *fname, const char *name, int m, int n, const double *A, int lda);

extern LIBSBBLAS_EXPORT void print_sdense2mm(FILE *f, const char *name, int m, int n, const float *A, int lda);
extern LIBSBBLAS_EXPORT void print_ddense2mm(FILE *f, const char *name, int m, int n, const double *A, int lda);

extern LIBSBBLAS_EXPORT void fprint_suff_sdense2mm(const char *fname, int suff, const char *name, int m, int n, const float *A, int lda);  
extern LIBSBBLAS_EXPORT void fprint_suff_ddense2mm(const char *fname, int suff, const char *name, int m, int n, const double *A, int lda);

extern LIBSBBLAS_EXPORT void print_shb(FILE *f, const char *name, hbmat_t *A, int full);
extern LIBSBBLAS_EXPORT void print_dhb(FILE *f, const char *name, hbmat_t *A, int full);


#endif // __MATFPRINT_H__
