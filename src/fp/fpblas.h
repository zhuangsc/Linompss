#ifndef __FPBLAS_H__
#define __FPBLAS_H__


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#define max(a,ts) ((a) >= (ts) ? (a) : (ts))
#define min(a,ts) ((a) >= (ts) ? (ts) : (a))


#if USE_MKL

#include "mkl.h"
#include "fpmatr.h"


#ifdef SINGLE_PRECISION

#define BLAS_cp(n, dx, incx, dy, incy) 												cblas_scopy(n, dx, incx, dy, incy)
#define BLAS_copy(n, x, incx, y, incy)												cblas_scopy(n, x, incx, y, incy)
#define BLAS_gemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) 	cblas_sgemm(CblasColMajor, transa.e, transb.e, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
#define BLAS_dot(n, dx, incx, dy, incy) 											cblas_sdot(n, dx, incx, dy, incy) 
#define BLAS_axpy(n, da, dx, incx, dy, incy) 										cblas_saxpy(n, da, dx, incx, dy, incy)
#define BLAS_scal(n, da, dx, incx)													cblas_sscal(n, da, dx, incx)
#define BLAS_nrm2(n, x, incx)														cblas_snrm2(n, x, incx)
#define BLAS_gemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)           	cblas_sgemv(CblasColMajor, trans.e, m, n, alpha, A, lda, x, incx, beta, y, incy)
#define BLAS_syrk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc) 					cblas_ssyrk(CblasColMajor, uplo.e, trans.e, n, k, alpha, A, lda, beta, C, ldc)
#define BLAS_trsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb) 			cblas_strsm(CblasColMajor, side.e, uplo.e, transa.e, diag.e, m, n, alpha, A, lda, B, ldb)
#define BLAS_trmm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb) 			cblas_strmm(CblasColMajor, side.e, uplo.e, transa.e, diag.e, m, n, alpha, A, lda, B, ldb)

#endif

#ifdef DOUBLE_PRECISION

#define BLAS_cp(n, dx, incx, dy, incy) 												cblas_dcopy(n, dx, incx, dy, incy)
#define BLAS_copy(n, x, incx, y, incy)												cblas_dcopy(n, x, incx, y, incy)
#define BLAS_gemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) 	cblas_dgemm(CblasColMajor, transa.e, transb.e, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
#define BLAS_dot(n, dx, incx, dy, incy) 											cblas_ddot(n, dx, incx, dy, incy) 
#define BLAS_axpy(n, da, dx, incx, dy, incy) 										cblas_daxpy(n, da, dx, incx, dy, incy)
#define BLAS_scal(n, da, dx, incx)													cblas_dscal(n, da, dx, incx)
#define BLAS_nrm2(n, x, incx)														cblas_dnrm2(n, x, incx)
#define BLAS_gemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)           	cblas_dgemv(CblasColMajor, trans.e, m, n, alpha, A, lda, x, incx, beta, y, incy)
#define BLAS_syrk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc) 					cblas_dsyrk(CblasColMajor, uplo.e, trans.e, n, k, alpha, A, lda, beta, C, ldc)
#define BLAS_trsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb) 			cblas_dtrsm(CblasColMajor, side.e, uplo.e, transa.e, diag.e, m, n, alpha, A, lda, B, ldb)
#define BLAS_trmm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb) 			cblas_dtrmm(CblasColMajor, side.e, uplo.e, transa.e, diag.e, m, n, alpha, A, lda, B, ldb)

#endif 


#else /* NOT MKL */


#include "blas.h"
#include "fpmatr.h"


#ifdef SINGLE_PRECISION

#define BLAS_cp(n, dx, incx, dy, incy) 												scopy_(&n, dx, &incx, dy, &incy)
#define BLAS_gemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) 	sgemm_(transa.s, transb.s, &(m), &(n), &(k), &(alpha), A, &(lda), B, &(ldb), &(beta), C, &(ldc))
#define BLAS_dot(n, dx, incx, dy, incy) 											sdot_(&(n), dx, &(incx), dy, &(incy)) 
#define BLAS_axpy(n, da, dx, incx, dy, incy) 										saxpy_(&n, &da, dx, &incx, dy, &incy)
#define BLAS_scal(n, da, dx, incx)													sscal_(&n, &da, dx, &incx)
#define BLAS_nrm2(n, x, incx)														snrm2_(&n, x, &incx)
#define BLAS_copy(n, x, incx, y, incy)												scopy_(&n, x, &incx, y, &incy)
#define BLAS_gemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy) 				sgemv_(trans.s, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy)
#define BLAS_trsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb) 			strsm_(side.s, uplo.s, transa.s, diag.s, &m, &n, &alpha, A, &lda, B, &ldb)
#define BLAS_syrk(trian, trans, n, k, alpha, A, lda, beta, C, ldc)					ssyrk_(trian.s, trans.s, &n, &k, &alpha, A, &lda, &beta, C, &ldc)
#define BLAS_trmm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb) 			strmm_(side.s, uplo.s, transa.s, diag.s, &m, &n, &alpha, A, &lda, B, &ldb)

#else

#define BLAS_cp(n, dx, incx, dy, incy) 												dcopy_(&(n), dx, &(incx), dy, &(incy))
#define BLAS_gemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) 	dgemm_(transa.s, transb.s, &(m), &(n), &(k), &(alpha), A, &(lda), B, &(ldb), &(beta), C, &(ldc))
#define BLAS_dot(n, dx, incx, dy, incy) 											ddot_(&(n), dx, &(incx), dy, &(incy)) 
#define BLAS_axpy(n, da, dx, incx, dy, incy) 										daxpy_(&n, &da, dx, &incx, dy, &incy)
#define BLAS_scal(n, da, dx, incx)													dscal_(&n, &da, dx, &incx)
#define BLAS_nrm2(n, x, incx)														dnrm2_(&n, x, &incx)
#define BLAS_copy(n, x, incx, y, incy)												dcopy_(&n, x, &incx, y, &incy)
#define BLAS_gemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy) 				dgemv_(trans.s, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy)
#define BLAS_trsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb) 			dtrsm_(side.s, uplo.s, transa.s, diag.s, &m, &n, &alpha, A, &lda, B, &ldb)
#define BLAS_syrk(trian, trans, n, k, alpha, A, lda, beta, C, ldc)					dsyrk_(trian.s, trans.s, &n, &k, &alpha, A, &lda, &beta, C, &ldc)
#define BLAS_trmm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb) 			dtrmm_(side.s, uplo.s, transa.s, diag.s, &m, &n, &alpha, A, &lda, B, &ldb)


#endif // PRECISION


#endif // MKL


#endif // __FPBLAS_H__
