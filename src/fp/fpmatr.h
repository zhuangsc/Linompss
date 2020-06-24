#ifndef __FPMATR_H__
#define __FPMATR_H__


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#define UPPER			1
#define LOWER			0


typedef union { 
	int e;
	char *s;
} ompssblas_t;


extern ompssblas_t OMPSSBLAS_ROWMAJOR;
extern ompssblas_t OMPSSBLAS_COLMAJOR;
extern ompssblas_t OMPSSBLAS_TRANSP;
extern ompssblas_t OMPSSBLAS_NTRANSP;
extern ompssblas_t OMPSSBLAS_UPPERTRIANG;
extern ompssblas_t OMPSSBLAS_LOWERTRIANG;
extern ompssblas_t OMPSSBLAS_NTRIANG;
extern ompssblas_t OMPSSBLAS_DIAGUNIT;
extern ompssblas_t OMPSSBLAS_NDIAGUNIT;
extern ompssblas_t OMPSSBLAS_LEFT;
extern ompssblas_t OMPSSBLAS_RIGHT;

typedef union { 
	char e;
	char *s;
} ompsslapack_t;

extern ompsslapack_t OMPSSLAPACK_UPPERTRIANG;
extern ompsslapack_t OMPSSLAPACK_LOWERTRIANG;


#if USE_MKL
#include "mkl.h"
#else
#include "blas.h"
#endif

static inline __attribute__((always_inline)) int TEST_LEFT(ompssblas_t descr) 
{
#if USE_MKL
	return descr.e == CblasLeft;	
#else
	return descr.s[0] == 'L';
#endif
}

static inline __attribute__((always_inline)) int TEST_RIGHT(ompssblas_t descr) 
{
#if USE_MKL
	return descr.e == CblasRight;	
#else
	return descr.s[0] == 'R';
#endif
}

static inline __attribute__((always_inline)) int TEST_UPPER(ompssblas_t descr) 
{
#if USE_MKL
	return descr.e == CblasUpper;	
#else
	return descr.s[0] == 'U';
#endif
}

static inline __attribute__((always_inline)) int TEST_LOWER(ompssblas_t descr) 
{
#if USE_MKL
	return descr.e == CblasLower;	
#else
	return descr.s[0] == 'L';
#endif
}

static inline __attribute__((always_inline)) int TEST_TRANSP(ompssblas_t descr) 
{
#if USE_MKL
	return descr.e == CblasTrans;	
#else
	return descr.s[0] == 'T';
#endif
}

static inline __attribute__((always_inline)) int TEST_NTRANSP(ompssblas_t descr) 
{
#if USE_MKL
	return descr.e != CblasTrans;	
#else
	return descr.s[0] == 'N';
#endif
}

static inline __attribute__((always_inline)) int TEST_NDIAGUNIT(ompssblas_t descr) 
{
#if USE_MKL
	return descr.e == CblasNonUnit;	
#else
	return descr.s[0] == 'N';
#endif
}

static inline __attribute__((always_inline)) int TEST_DIAGUNIT(ompssblas_t descr) 
{
#if USE_MKL
	return descr.e == CblasUnit;	
#else
	return descr.s[0] == 'U';
#endif
}


#endif // __FPMATR_H__
