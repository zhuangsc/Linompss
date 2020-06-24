#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "task_axpy.h"

#include "fptype.h"
#include "fpblas.h"
#include "blas.h"
#include "blas_ext.h"
#include "selfsched.h"




#ifdef DOUBLE_PRECISION

#define __t_axpy				task_daxpy
#define __t_scal_axpy			task_scal_daxpy
#define __t_ext_axpy			task_ext_daxpy
#define __t_extm_axpy			task_extm_daxpy
#define __t_cpaxpy				task_dcpaxpy
#define __t_scal_cpaxpy			task_scal_dcpaxpy
#define __t_cpaxpy_comb			task_dcpaxpy_comb
#define __t_cpaxpy_comb4		task_dcpaxpy_comb4
#define	__t_axpy4				task_daxpy4

#else

#define __t_axpy				task_saxpy
#define __t_scal_axpy			task_scal_saxpy
#define __t_ext_axpy			task_ext_saxpy
#define __t_extm_axpy			task_extm_saxpy
#define __t_cpaxpy				task_scpaxpy
#define __t_scal_cpaxpy			task_scal_scpaxpy
#define __t_cpaxpy_comb			task_scpaxpy_comb
#define __t_cpaxpy_comb4		task_scpaxpy_comb4
#define	__t_axpy4				task_saxpy4

#endif


void __t_axpy(int bm, int bn, int m, int n, fp_t *Anum, fp_t *Aden, fp_t *X, fp_t *Y) {
#ifndef OMPSS_DRYRUN
	int j;
	for ( j=0; j<bn; ++j) {
		fp_t f = Anum[j] / Aden[j];
		BLAS_axpy(bm, f, X, i_one, Y, i_one);
		X += m;
		Y += m;
	}
#endif
}


void __t_scal_axpy(int bm, int bn, int m, int n, fp_t alpha, fp_t *Anum, fp_t *Aden, fp_t *X, fp_t *Y) {
#ifndef OMPSS_DRYRUN
	fp_t localA[bn];
	int j;
	for ( j = 0; j < bn; ++j ) {	
		localA[j] = alpha * Anum[j] / Aden[j];
	}

	for ( j=0; j<bn; ++j) {
		BLAS_axpy(bm, localA[j], X, i_one, Y, i_one);
		X += m;
		Y += m;
	}
#endif
}


void __t_ext_axpy(int bm, int bn, int m, int n, fp_t SA, fp_t *X, fp_t *Y, fp_t *Z) {
#ifndef OMPSS_DRYRUN
	//verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, X);
	//verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, Y);
	int j;
	for ( j=0; j<bn; ++j) {
		BLAS_cp(bm, &Y[j*m], i_one, &Z[j*m], i_one);
		BLAS_axpy(bm, SA, &X[j*m], i_one, &Z[j*m], i_one);
	}
	//verify_2df_array(__FILE__, __LINE__, M, N, BSA, BSB, Z);
#endif
}


void __t_extm_axpy(int bm, int bn, int m, int n, fp_t *SAnum, fp_t *SAden, fp_t *X, fp_t *Y, fp_t *Z, int p) 
{
	int j;
	for ( j=0; j<bn; ++j) {
		fp_t f = SAnum[j] / SAden[j];
		BLAS_cp(bm, &Y[j*m], i_one, &Z[j*m], i_one);
		BLAS_axpy(bm, f, &X[j*m], i_one, &Z[j*m], i_one);
	}
}


#ifdef DOUBLE_PRECISION
void task_sdaxpy(int bm, int bn, int m, int n, float a, float *x, double *y) {
	int i;
	for ( i=0; i<bm; ++i ) {
		int j;
		for ( j=0; j<bn; ++j ) {
			y[j*m+i] += a * x[j*m+i];
		}
	}
}

void task_sdcpaxpy(int bm, int bn, int m, int n, double a, double *x, double *y, double *z) {
	int i;
	for ( i=0; i<bm; ++i ) {
		int j;
		for ( j=0; j<bn; ++j ) {
			z[j*m+i] = y[j*m+i] + a * x[j*m+i];
		}
	}
}

#endif


void __t_cpaxpy(int bm, int bn, int m, int n, fp_t *Anum, fp_t *Aden, fp_t *X, fp_t *Y, fp_t *Z) {
	int j;
	for ( j=0; j<bn; ++j) {
		fp_t factor = Anum[j] / Aden[j];
		BLAS_cp(bm, Y, i_one, Z, i_one);
		BLAS_axpy(bm, factor, X, i_one, Z, i_one);
		X += m;
		Y += m;
		Z += m;
	}
}


void __t_scal_cpaxpy(int bm, int bn, int m, int n, fp_t alpha, fp_t *Anum, fp_t *Aden, fp_t *X, fp_t *Y, fp_t *Z) {
	fp_t localA[bn];
	int j;
	for ( j = 0; j < bn; ++j ) {	
		localA[j] = alpha * Anum[j] / Aden[j];
	}

	for ( j=0; j<bn; ++j) {
		BLAS_cp(bm, Y, i_one, Z, i_one);
		BLAS_axpy(bm, localA[j], X, i_one, Z, i_one);
		X += m;
		Y += m;
		Z += m;
	}
}


void __t_cpaxpy_comb(int bm, int bn, int m, int n, fp_t alpha, fp_t *Anum, fp_t *Aden, fp_t *X1, fp_t *X2, fp_t *Y1, fp_t *Y2, fp_t *Z1, fp_t *Z2) 
{
	int j;
	for ( j=0; j<bn; ++j) {
		/* update of x */
		fp_t factor = Anum[j] / Aden[j];
//		printf("alpha1: %e, alpha2: %e, factor: %e\n", Anum[j], Aden[j], factor);
		BLAS_cp(bm, Y2, i_one, Z2, i_one);
		BLAS_axpy(bm, factor, X2, i_one, Z2, i_one);
		X2 += m;
		Y2 += m;
		Z2 += m;

		
		/* update of r */
		factor = alpha * factor;
		BLAS_cp(bm, Y1, i_one, Z1, i_one);
		BLAS_axpy(bm, factor, X1, i_one, Z1, i_one);
		X1 += m;
		Y1 += m;
		Z1 += m;
	}
}


void __t_cpaxpy_comb4(int bm, int bn, int m, int n, fp_t *gamma1, fp_t *gamma2, fp_t *delta, fp_t *sigma2, fp_t *P2, fp_t *V2, fp_t *X2, fp_t *R2, fp_t *S,\
	fp_t *sigma1, fp_t *P1, fp_t *V1, fp_t *X1, fp_t *R1) {
	int j;
	for ( j=0; j<bn; ++j) {
		fp_t beta = gamma1[j] / gamma2[j];
		fp_t lsigma1 = sigma1[j] = delta[j] - beta * beta * sigma2[j];

		/* update of P1 */
		BLAS_cp(bm, R2, i_one, P1, i_one);
		BLAS_axpy(bm, beta, P2, i_one, P1, i_one);

		/* update of V1 */
		BLAS_cp(bm, S, i_one, V1, i_one);
		BLAS_axpy(bm, beta, V2, i_one, V1, i_one);
		
		fp_t alpha = gamma1[j] / lsigma1;
		//printf("cpaxpy_comb4: gamma %.16e %.16e sigma %.16e beta %.16e alpha %.16e\n", gamma1[j], gamma2[j], lsigma1, beta, alpha);
		/* update of X1 */
		BLAS_cp(bm, X2, i_one, X1, i_one);
		BLAS_axpy(bm, alpha, P1, i_one, X1, i_one);

		/* update of R1 */
		alpha = - alpha;
		BLAS_cp(bm, R2, i_one, R1, i_one);
		BLAS_axpy(bm, alpha, V1, i_one, R1, i_one);

		P2 += m;
		V2 += m;
		X2 += m;
		R2 += m;
		S += m;
		P1 += m;
		V1 += m;
		X1 += m;
		R1 += m;
	}
}


void __t_axpy4(int bm, int bn, int m, int n, fp_t *gamman, fp_t *gammad, fp_t *alpha, fp_t *z, fp_t *s, fp_t *p, fp_t *q, fp_t *w, fp_t *xp, fp_t *x, fp_t *rp, fp_t *r) 
{
	fp_t *tmp = malloc(bm * sizeof(fp_t));

	int j;
	for ( j=0; j<bn; ++j) {
		int init = gammad[j] < FP_NOUGHT;
		fp_t beta = init ? FP_NOUGHT: gamman[j] / gammad[j];
		fp_t lalpha = alpha[j];
		fp_t lmalpha = -lalpha;

		//printf("beta %f gammap %f\n", beta, gammad[j]);
		/* z = q + beta * z */
		BLAS_cp(bm, z, I_ONE, tmp, I_ONE);
		/* warning! */
		BLAS_set(FP_NOUGHT, bm, 1, z, m);
		BLAS_axpy(bm, beta, tmp, I_ONE, z, I_ONE);
		BLAS_axpy(bm, FP_ONE, q, I_ONE, z, I_ONE);

		/* s = w + beta * s */
		BLAS_cp(bm, s, I_ONE, tmp, I_ONE);
		BLAS_set(FP_NOUGHT, bm, 1, s, m);
		BLAS_axpy(bm, beta, tmp, I_ONE, s, I_ONE);
		BLAS_axpy(bm, FP_ONE, w, I_ONE, s, I_ONE);

		/* p = r + beta * p */
		BLAS_cp(bm, p, I_ONE, tmp, I_ONE);
		BLAS_set(FP_NOUGHT, bm, 1, p, m);
		BLAS_axpy(bm, beta, tmp, I_ONE, p, I_ONE);
		BLAS_axpy(bm, FP_ONE, r, I_ONE, p, I_ONE);

		/* x = x + alpha * p */
		BLAS_cp(bm, xp, I_ONE, x, I_ONE);
		BLAS_axpy(bm, lalpha, p, I_ONE, x, I_ONE);
		/* r = r - alpha * s */
		BLAS_cp(bm, rp, I_ONE, r, I_ONE);
		BLAS_axpy(bm, lmalpha, s, I_ONE, r, I_ONE);
		/* w = w - alpha * z */
		BLAS_axpy(bm, lmalpha, z, I_ONE, w, I_ONE);

		z += m;
		s += m;
		p += m;
		q += m;
		w += m;
		xp += m;
		x += m;
		rp += m;
		r += m;
	}

	free(tmp);
}

