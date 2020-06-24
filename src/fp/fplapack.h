#ifndef __FPLAPACK_H__
#define __FPLAPACK_H__


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "fpblas.h"


#if USE_MKL

#include "mkl.h"

typedef char fplapackd_t;

#define LAPACK_TRIANG_UPPER								'U'
#define LAPACK_TRIANG_LOWER								'L'
#define LAPACK_LEFT										'L'
#define LAPACK_RIGHT									'R'
#define LAPACK_TRANSP									'T'
#define LAPACK_NOTRANSP									'N'
#define LAPACK_FORWARD									'F'
#define LAPACK_BACKWARD									'B'
#define LAPACK_COLWISE									'C'
#define LAPACK_ROWWISE									'R'


#ifdef SINGLE_PRECISION

#define LAPACK_getrf(m, n, A, lda, ipiv, info) 												LAPACKE_sgetrf(LAPACK_COL_MAJOR, m, n, A, lda, ipiv )
#define LAPACK_laswp(n, A, lda, k1, k2, ipiv, incx)											LAPACKE_slaswp(LAPACK_COL_MAJOR, n, A, lda, k1, k2, ipiv, incx)
#define LAPACK_potrf(trian, n, A, lda, info)												LAPACKE_spotrf(LAPACK_COL_MAJOR, trian.e, n, A, lda)
#define LAPACK_geqrt(m, n, b, A, lda, T, ldt, work, info)									LAPACKE_sgeqrt(LAPACK_COL_MAJOR, m, n, b, A, lda, T, ldt)
#define LAPACK_lacpy(uplo, m, n, A, lda, B, ldb)											LAPACKE_slacpy(LAPACK_COL_MAJOR, uplo, m, n, A, lda, B, ldb)
// MKL in MareNostrum does not have a correct declaration in header file. 
#define LAPACK_larfb(side, trans, direct, storev, m, n, k, V, ldv, T, ldt, C, ldc, W, ldw)	LAPACKE_slarfb(LAPACK_COL_MAJOR, side, trans, direct, storev, m, n, k, V, ldv, T, ldt, C, ldc)
#define LAPACK_larft(direct, storev, n, k, V, ldv, tau, T, ldt)								LAPACKE_slarft(LAPACK_COL_MAJOR, direct, storev, n, k, V, ldv, tau, T, ldt)
#define LAPACK_geqr2(m, n, A, lda, tau, W, info )											LAPACKE_sgeqr2(LAPACK_COL_MAJOR, m, n, A, lda, tau)//, W, info )
#define LAPACK_lange(norm, m, n, A, lda)													LAPACKE_slange(LAPACK_COL_MAJOR, norm, m, n, A, lda)

#else

#define LAPACK_getrf(m, n, A, lda, ipiv, info) 												LAPACKE_dgetrf(LAPACK_COL_MAJOR, m, n, A, lda, ipiv )
#define LAPACK_laswp(n, A, lda, k1, k2, ipiv, incx)											LAPACKE_dlaswp(LAPACK_COL_MAJOR, n, A, lda, k1, k2, ipiv, incx)
#define LAPACK_potrf(trian, n, A, lda, info)												LAPACKE_dpotrf(LAPACK_COL_MAJOR, trian.e, n, A, lda)
#define LAPACK_geqrt(m, n, b, A, lda, T, ldt, work, info)									LAPACKE_dgeqrt(LAPACK_COL_MAJOR, m, n, b, A, lda, T, ldt)
#define LAPACK_lacpy(uplo, m, n, A, lda, B, ldb)											LAPACKE_dlacpy(LAPACK_COL_MAJOR, uplo, m, n, A, lda, B, ldb)
#define LAPACK_larfb(side, trans, direct, storev, m, n, k, V, ldv, T, ldt, C, ldc, W, ldw)	LAPACKE_dlarfb(LAPACK_COL_MAJOR, side, trans, direct, storev, m, n, k, V, ldv, T, ldt, C, ldc)
#define LAPACK_larft(direct, storev, n, k, V, ldv, tau, T, ldt)								LAPACKE_dlarft(LAPACK_COL_MAJOR, direct, storev, n, k, V, ldv, tau, T, ldt)
#define LAPACK_geqr2(m, n, A, lda, tau, W, info )											LAPACKE_dgeqr2(LAPACK_COL_MAJOR, m, n, A, lda, tau)//, W, info )
#define LAPACK_lange(norm, m, n, A, lda)													LAPACKE_dlange(LAPACK_COL_MAJOR, norm, m, n, A, lda)

#endif 


#else /* NOT MKL */


#include "blas.h"

typedef const char * fplapackd_t;

#define LAPACK_TRIANG_UPPER								"U"
#define LAPACK_TRIANG_LOWER								"L"
#define LAPACK_LEFT										"L"
#define LAPACK_RIGHT									"R"
#define LAPACK_TRANSP									"T"
#define LAPACK_NOTRANSP									"N"
#define LAPACK_FORWARD									"F"
#define LAPACK_BACKWARD									"B"
#define LAPACK_COLWISE									"C"
#define LAPACK_ROWWISE									"R"


#ifdef SINGLE_PRECISION

#define LAPACK_getrf(m, n, A, lda, ipiv, info) 												sgetrf_(&m, &n, A, &lda, ipiv, info)
#define LAPACK_laswp(n, A, lda, k1, k2, ipiv, incx)											slaswp_(&n, A, &lda, &k1, &k2, ipiv, &incx)
#define LAPACK_potrf(trian, n, A, lda, info)												spotrf_(trian.s, &n, A, &lda, &info)
#define LAPACK_geqrt(m, n, b, A, lda, T, ldt, work, info)									sgeqrt_(&m, &n, &b, A, &lda, T, &ldt, work, &info)
#define LAPACK_lacpy(uplo, m, n, A, lda, B, ldb)											slacpy_(uplo, &m, &n, A, &lda, B, &ldb)
#define LAPACK_larfb(side, trans, direct, storev, m, n, k, V, ldv, T, ldt, C, ldc, W, ldw)	slarfb_(side, trans, direct, storev, &m, &n, &k, V, &ldv, T, &ldt, C, &ldc, W, &ldw)
#define LAPACK_geqr2(m, n, A, lda, tau, W, info )											sgeqr2_(&m, &n, A, &lda, tau, W, &info )
#define LAPACK_larft(direct, storev, n, k, V, ldv, tau, T, ldt)								slarft_(direct, storev, &n, &k, V, &ldv, tau, T, &ldt)
#warning Zhuang, those *lange functions expect a work array
#define LAPACK_lange(norm, m, n, A, lda)													slange_(norm, m, n, A, lda, NULL)

#else

#define LAPACK_getrf(m, n, A, lda, ipiv, info) 												dgetrf_(&m, &n, A, &lda, ipiv, info)
#define LAPACK_laswp(n, A, lda, k1, k2, ipiv, incx)											dlaswp_(&n, A, &lda, &k1, &k2, ipiv, &incx)
#define LAPACK_potrf(trian, n, A, lda, info)												dpotrf_(trian.s, &n, A, &lda, &info)
#define LAPACK_geqrt(m, n, b, A, lda, T, ldt, work, info)									dgeqrt_(&m, &n, &b, A, &lda, T, &ldt, work, &info)
#define LAPACK_lacpy(uplo, m, n, A, lda, B, ldb)											dlacpy_(uplo, &m, &n, A, &lda, B, &ldb)
#define LAPACK_larfb(side, trans, direct, storev, m, n, k, V, ldv, T, ldt, C, ldc, W, ldw)	dlarfb_(side, trans, direct, storev, &m, &n, &k, V, &ldv, T, &ldt, C, &ldc, W, &ldw)
#define LAPACK_geqr2(m, n, A, lda, tau, W, info )											dgeqr2_(&m, &n, A, &lda, tau, W, &info )
#define LAPACK_larft(direct, storev, n, k, V, ldv, tau, T, ldt)								dlarft_(direct, storev, &n, &k, V, &ldv, tau, T, &ldt)
#warning Zhuang, those *lange functions expect a work array
#define LAPACK_lange(norm, m, n, A, lda)													dlange_(norm, m, n, A, lda, NULL)

#endif // PRECISION


#endif // MKL


#endif // __FPLAPACK_H__
