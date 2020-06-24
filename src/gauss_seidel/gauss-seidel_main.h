#ifndef __JACOBI_MAIN_H__
#define __JACOBI_MAIN_H__


#include "fptype.h"
#include "hb.h"


#ifdef SINGLE_PRECISION

#define GS_MAIN_CSR		sgs_main_csr

#else

#define GS_MAIN_CSR		dgs_main_csr

#endif

int sgs_main_csr(hbmat_t *Acsr, float *x, float *b, int bs, int max_iter, hbmat_t **diagL, int lookahead, float threshold, float *work);

int dgs_main_csr(hbmat_t *Acsr, double *x, double *b, int bs, int max_iter, hbmat_t **diagL, int lookahead, double threshold, double *work);

#endif
