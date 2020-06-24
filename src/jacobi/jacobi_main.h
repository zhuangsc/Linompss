#ifndef __JACOBI_MAIN_H__
#define __JACOBI_MAIN_H__


#include "fptype.h"
#include "hb.h"

#ifdef SINGLE_PRECISION
#define JACOBI_MAIN_CSR		sjacobi_main_csr
#else
#define JACOBI_MAIN_CSR		djacobi_main_csr
#endif

//unsigned int jacobi_main_csr(hbmat_t *Acsr, fp_t *x, fp_t *b, int bs, int max_iter, hbmat_t **diagL, int lookahead, fp_t res, fp_t *work);

unsigned int sjacobi_main_csr(hbmat_t *Acsr, float *x, float *b, int bs, int max_iter, hbmat_t **diagL, int lookahead, float res, float *work, int *res_p);

unsigned int djacobi_main_csr(hbmat_t *Acsr, double *x, double *b, int bs, int max_iter, hbmat_t **diagL, int lookahead, double res, double *work, int *res_p);

#endif
