#ifndef __FPSBLAS_H__
#define __FPSBLAS_H__


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#if USE_MKL

#include "mkl.h"
#include "fpmatr.h"

typedef int fpblasd_t;

#ifdef SINGLE_PRECISION

#define SBLAS_csrmv(trans, m, n, alpha, matdescra, vval, vpos ,vptr, vptr1, X, beta, Y)  			mkl_scsrmv(trans, &m, &n, &alpha, matdescra, vval, vpos ,vptr, vptr1, X, &beta, Y)
#define SBLAS_csrsv(trans, m, alpha, matdescra, vval, vpos, vptr, vptr1, X, Y) 						mkl_scsrsv(trans, m, alpha, matdescra, vval, vpos, vptr, vptr1, X, Y)
#define SBLAS_CSRCSC(job, m, acsr, ja, ia, acsc, ja1, ia1, info)									mkl_scsrcsc(job, m, acsr, ja, ia, acsc, ja1, ia1, info)
#define SBLAS_csrmm(trans, m, n, k, alpha, matdescra, val, vpos, vptr, vptr1, B, ldb, beta, C, ldc)	mkl_scsrmm(trans, &m, &n, &k, &alpha, matdescra, val, vpos, vptr, vptr1, B, &ldb, &beta, C, &ldc)
#define SBLAS_csrsm(trans, m, n, alpha, matdescra, vval, vpos, vptr, vptr1, B, ldb, C, ldc)			mkl_scsrsm(trans, &m, &n, &alpha, matdescra, vval, vpos, vptr, vptr1, B, &ldb, C, &ldc)

#else

#define SBLAS_csrmv(trans, m, n, alpha, matdescra, vval, vpos ,vptr, vptr1, X, beta, Y)  			mkl_dcsrmv(trans, &m, &n, &alpha, matdescra, vval, vpos ,vptr, vptr1, X, &beta, Y)
#define SBLAS_csrsv(trans, m, alpha, matdescra, vval, vpos, vptr, vptr1, X, Y) 						mkl_dcsrsv(trans, m, alpha, matdescra, vval, vpos, vptr, vptr1, X, Y)
#define SBLAS_CSRCSC(job, m, acsr, ja, ia, acsc, ja1, ia1, info)									mkl_dcsrcsc(job, m, acsr, ja, ia, acsc, ja1, ia1, info)
#define SBLAS_csrmm(trans, m, n, k, alpha, matdescra, val, vpos, vptr, vptr1, B, ldb, beta, C, ldc)	mkl_dcsrmm(trans, &m, &n, &k, &alpha, matdescra, val, vpos, vptr, vptr1, B, &ldb, &beta, C, &ldc)
#define SBLAS_csrsm(trans, m, n, alpha, matdescra, vval, vpos, vptr, vptr1, B, ldb, C, ldc)			mkl_dcsrsm(trans, &m, &n, &alpha, matdescra, vval, vpos, vptr, vptr1, B, &ldb, C, &ldc)

#endif 


#else /* NOT MKL */


#ifdef SINGLE_PRECISION

#define SBLAS_csrmv(trans, m, n, alpha, matdescra, vval, vpos ,vptr, vptr1, X, beta, Y)  			
#define SBLAS_csrsv(trans, m, alpha, matdescra, vval, vpos, vptr, vptr1, X, Y) 						
#define SBLAS_CSRCSC(job, m, acsr, ja, ia, acsc, ja1, ia1, info)									
#define SBLAS_csrmm(trans, m, n, k, alpha, matdescra, val, vpos, vptr, vptr1, B, ldb, beta, C, ldc)	

#else

#define SBLAS_csrmv(trans, m, n, alpha, matdescra, vval, vpos ,vptr, vptr1, X, beta, Y)  			
#define SBLAS_csrsv(trans, m, alpha, matdescra, vval, vpos, vptr, vptr1, X, Y) 						
#define SBLAS_CSRCSC(job, m, acsr, ja, ia, acsc, ja1, ia1, info)									
#define SBLAS_csrmm(trans, m, n, k, alpha, matdescra, val, vpos, vptr, vptr1, B, ldb, beta, C, ldc)	

#endif 


#endif // MKL


#endif // __FPSBLAS_H__
