#ifndef __BLAS_H__
#define __BLAS_H__


#define F77_CALL(x)  x ## _
#define F77_NAME(x)  F77_CALL(x)

#ifdef  __cplusplus
extern "C" {
#endif

/* Level 1 BLAS */

extern double F77_NAME(dasum)(const int *n, const double *dx, const int *incx);
extern void   F77_NAME(daxpy)(const int *n, const double *alpha,
                              const double *dx, const int *incx,
                              double *dy, const int *incy);
extern void   F77_NAME(dcopy)(const int *n, const double *dx, const int *incx,
                              double *dy, const int *incy);
extern double F77_NAME(ddot) (const int *n, const double *dx, const int *incx,
                              const double *dy, const int *incy);
extern double F77_NAME(dnrm2)(const int *n, const double *dx, const int *incx);
extern void   F77_NAME(drot) (const int *n, double *dx, const int *incx,
                              double *dy, const int *incy, const double *c,
                              const double *s);
extern void   F77_NAME(drotg)(const double *a, const double *b, double *c, 
                              double *s);
extern void   F77_NAME(drotm)(const int *n, double *dx, const int *incx,
                              double *dy, const int *incy,const double *dparam);
extern void   F77_NAME(drotmg)(const double *dd1, const double *dd2, 
                               const double *dx1, const double *dy1, 
                               double *param);
extern void   F77_NAME(dscal)(const int *n, const double *alpha, double *dx,
                              const int *incx);
extern void   F77_NAME(dswap)(const int *n, double *dx, const int *incx,
                              double *dy, const int *incy);
extern int    F77_NAME(idamax)(const int *n, const double *dx, const int *incx);


extern float F77_NAME(sasum)(const int *n, const float *dx, const int *incx);
extern void   F77_NAME(saxpy)(const int *n, const float *alpha,
                              const float *dx, const int *incx,
                              float *dy, const int *incy);
extern void   F77_NAME(scopy)(const int *n, const float *dx, const int *incx,
                              float *dy, const int *incy);
extern float F77_NAME(sdot) (const int *n, const float *dx, const int *incx,
                              const float *dy, const int *incy);
extern float F77_NAME(snrm2)(const int *n, const float *dx, const int *incx);
extern void   F77_NAME(srot) (const int *n, float *dx, const int *incx,
                              float *dy, const int *incy, const float *c,
                              const float *s);
extern void   F77_NAME(srotg)(const float *a, const float *b, float *c, 
                              float *s);
extern void   F77_NAME(srotm)(const int *n, float *dx, const int *incx,
                              float *dy, const int *incy,const float *dparam);
extern void   F77_NAME(srotmg)(const float *dd1, const float *dd2, 
                               const float *dx1, const float *dy1, 
                               float *param);
extern void   F77_NAME(sscal)(const int *n, const float *alpha, float *dx,
                              const int *incx);
extern void   F77_NAME(sswap)(const int *n, float *dx, const int *incx,
                              float *dy, const int *incy);
extern int    F77_NAME(isamax)(const int *n, const float *dx, const int *incx);



/* Level 2 BLAS */

extern void   F77_NAME(dgbmv)(const char *trans, const int *m, const int *n,
                              const int *kl,const int *ku, const double *alpha,
                              const double *a, const int *lda, const double *x,
                              const int *incx, const double *beta, double *y,
                              const int *incy);
extern void   F77_NAME(dgemv)(const char *trans, const int *m, const int *n,
                              const double *alpha, const double *a,
                              const int *lda, const double *x, const int *incx,
                              const double *beta, double *y, const int *incy);
extern void   F77_NAME(dsbmv)(const char *uplo, const int *n, const int *k,
                              const double *alpha, const double *a,
                              const int *lda, const double *x, const int *incx,
                              const double *beta, double *y, const int *incy);
extern void   F77_NAME(dspmv)(const char *uplo, const int *n,
                              const double *alpha, const double *ap,
                              const double *x, const int *incx,
                              const double *beta, double *y, const int *incy);
extern void   F77_NAME(dsymv)(const char *uplo, const int *n,
                              const double *alpha, const double *a,
                              const int *lda, const double *x, const int *incx,
                              const double *beta, double *y, const int *incy);
extern void   F77_NAME(dtbmv)(const char *uplo, const char *trans,
                              const char *diag, const int *n, const int *k,
                              const double *a, const int *lda,
                              double *x, const int *incx);
extern void   F77_NAME(dtpmv)(const char *uplo, const char *trans,
                              const char *diag, const int *n, const double *ap,
                              double *x, const int *incx);
extern void   F77_NAME(dtrmv)(const char *uplo, const char *trans,
                              const char *diag, const int *n, const double *a,
                              const int *lda, double *x, const int *incx);
extern void   F77_NAME(dtbsv)(const char *uplo, const char *trans,
                              const char *diag, const int *n, const int *k,
                              const double *a, const int *lda,
                              double *x, const int *incx);
extern void   F77_NAME(dtpsv)(const char *uplo, const char *trans,
                              const char *diag, const int *n,
                              const double *ap, double *x, const int *incx);
extern void   F77_NAME(dtrsv)(const char *uplo, const char *trans,
                              const char *diag, const int *n,
                              const double *a, const int *lda,
                              double *x, const int *incx);
extern void   F77_NAME(dger) (const int *m, const int *n, const double *alpha,
                              double *x, const int *incx,
                              double *y, const int *incy,
                              double *a, const int *lda);
extern void   F77_NAME(dsyr) (const char *uplo, const int *n,
                              const double *alpha, const double *x,
                              const int *incx, double *a, const int *lda);
extern void   F77_NAME(dspr) (const char *uplo, const int *n,
                              const double *alpha, const double *x,
                              const int *incx, double *ap);
extern void   F77_NAME(dsyr2)(const char *uplo, const int *n, 
                              const double *alpha, const double *x,
                              const int *incx, const double *y, const int *incy,
                              double *a, const int *lda);
extern void   F77_NAME(dspr2)(const char *uplo, const int *n,
                              const double *alpha, const double *x,
                              const int *incx, const double *y,
                              const int *incy, double *ap);


extern void   F77_NAME(sgbmv)(const char *trans, const int *m, const int *n,
                              const int *kl,const int *ku, const float *alpha,
                              const float *a, const int *lda, const float *x,
                              const int *incx, const float *beta, float *y,
                              const int *incy);
extern void   F77_NAME(sgemv)(const char *trans, const int *m, const int *n,
                              const float *alpha, const float *a,
                              const int *lda, const float *x, const int *incx,
                              const float *beta, float *y, const int *incy);
extern void   F77_NAME(ssbmv)(const char *uplo, const int *n, const int *k,
                              const float *alpha, const float *a,
                              const int *lda, const float *x, const int *incx,
                              const float *beta, float *y, const int *incy);
extern void   F77_NAME(sspmv)(const char *uplo, const int *n,
                              const float *alpha, const float *ap,
                              const float *x, const int *incx,
                              const float *beta, float *y, const int *incy);
extern void   F77_NAME(ssymv)(const char *uplo, const int *n,
                              const float *alpha, const float *a,
                              const int *lda, const float *x, const int *incx,
                              const float *beta, float *y, const int *incy);
extern void   F77_NAME(stbmv)(const char *uplo, const char *trans,
                              const char *diag, const int *n, const int *k,
                              const float *a, const int *lda,
                              float *x, const int *incx);
extern void   F77_NAME(stpmv)(const char *uplo, const char *trans,
                              const char *diag, const int *n, const float *ap,
                              float *x, const int *incx);
extern void   F77_NAME(strmv)(const char *uplo, const char *trans,
                              const char *diag, const int *n, const float *a,
                              const int *lda, float *x, const int *incx);
extern void   F77_NAME(stbsv)(const char *uplo, const char *trans,
                              const char *diag, const int *n, const int *k,
                              const float *a, const int *lda,
                              float *x, const int *incx);
extern void   F77_NAME(stpsv)(const char *uplo, const char *trans,
                              const char *diag, const int *n,
                              const float *ap, float *x, const int *incx);
extern void   F77_NAME(strsv)(const char *uplo, const char *trans,
                              const char *diag, const int *n,
                              const float *a, const int *lda,
                              float *x, const int *incx);
extern void   F77_NAME(sger) (const int *m, const int *n, const float *alpha,
                              float *x, const int *incx,
                              float *y, const int *incy,
                              float *a, const int *lda);
extern void   F77_NAME(ssyr) (const char *uplo, const int *n,
                              const float *alpha, const float *x,
                              const int *incx, float *a, const int *lda);
extern void   F77_NAME(sspr) (const char *uplo, const int *n,
                              const float *alpha, const float *x,
                              const int *incx, float *ap);
extern void   F77_NAME(ssyr2)(const char *uplo, const int *n, 
                              const float *alpha, const float *x,
                              const int *incx, const float *y, const int *incy,
                              float *a, const int *lda);
extern void   F77_NAME(sspr2)(const char *uplo, const int *n,
                              const float *alpha, const float *x,
                              const int *incx, const float *y,
                              const int *incy, float *ap);



/* Level 3 BLAS */

extern void   F77_NAME(dgemm)(const char *transa, const char *transb,
                              const int *m, const int *n, const int *k,
                              const double *alpha, const double *a,
                              const int *lda, const double *b, const int *ldb,
                              const double *beta, double *c, const int *ldc);
extern void   F77_NAME(dtrsm)(const char *side, const char *uplo,
                              const char *transa, const char *diag,
                              const int *m, const int *n, const double *alpha,
                              const double *a, const int *lda,
                              double *b, const int *ldb);
extern void   F77_NAME(dtrmm)(const char *side, const char *uplo,
                              const char *transa, const char *diag,
                              const int *m, const int *n, const double *alpha,
                              const double *a, const int *lda,
                              double *b, const int *ldb);
extern void   F77_NAME(dsymm)(const char *side, const char *uplo, const int *m,
                              const int *n, const double *alpha,
                              const double *a, const int *lda,
                              const double *b, const int *ldb,
                              const double *beta, double *c, const int *ldc);
extern void   F77_NAME(dsyrk)(const char *uplo, const char *trans,
                              const int *n, const int *k,
                              const double *alpha, const double *a,
                              const int *lda, const double *beta,
                              double *c, const int *ldc);
extern void   F77_NAME(dsyr2k)(const char *uplo, const char *trans,
                               const int *n, const int *k,
                               const double *alpha, const double *a,
                               const int *lda, const double *b, const int *ldb,
                               const double *beta, double *c, const int *ldc);


extern void   F77_NAME(sgemm)(const char *transa, const char *transb,
                              const int *m, const int *n, const int *k,
                              const float *alpha, const float *a,
                              const int *lda, const float *b, const int *ldb,
                              const float *beta, float *c, const int *ldc);
extern void   F77_NAME(strsm)(const char *side, const char *uplo,
                              const char *transa, const char *diag,
                              const int *m, const int *n, const float *alpha,
                              const float *a, const int *lda,
                              float *b, const int *ldb);
extern void   F77_NAME(strmm)(const char *side, const char *uplo,
                              const char *transa, const char *diag,
                              const int *m, const int *n, const float *alpha,
                              const float *a, const int *lda,
                              float *b, const int *ldb);
extern void   F77_NAME(ssymm)(const char *side, const char *uplo, const int *m,
                              const int *n, const float *alpha,
                              const float *a, const int *lda,
                              const float *b, const int *ldb,
                              const float *beta, float *c, const int *ldc);
extern void   F77_NAME(ssyrk)(const char *uplo, const char *trans,
                              const int *n, const int *k,
                              const float *alpha, const float *a,
                              const int *lda, const float *beta,
                              float *c, const int *ldc);
extern void   F77_NAME(ssyr2k)(const char *uplo, const char *trans,
                               const int *n, const int *k,
                               const float *alpha, const float *a,
                               const int *lda, const float *b, const int *ldb,
                               const float *beta, float *c, const int *ldc);


#ifdef  __cplusplus
}
#endif

#endif // __BLAS_H__ 
