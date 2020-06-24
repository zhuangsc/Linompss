#include "chols_warm.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <stdio.h>

#include "hb.h"
#include "fptype.h"
#if USE_MKL
#include "mkl.h"
#endif


#define DIM 256
#define FILL (DIM * 0.7)


void chol_warmup() {
#if USE_MKL
	srand48(time(0));
	hbmat_t *t = malloc(sizeof(hbmat_t));
	t->m = DIM; t->n = DIM;
	t->vdiag = NULL;
	int m = t->m;
	int alpha = 1; int beta = 1;
	int *vptr = t->vptr = malloc((DIM+1) * sizeof(int));
	int *vpos = t->vpos = malloc((DIM * DIM) *sizeof(int));
	fp_t *vval = t->vval = malloc((DIM*DIM)*sizeof(fp_t));
	vptr[0] = 0;
	int vpos_p = 0;
	puts("warm-up");
	for ( int i = 1; i <= DIM; ++i ) {
		vptr[i] = vptr[i-1] +  FILL;
		int vp = 0;
		for ( int j = vptr[i-1]; j < vptr[i]; ++j ) {
			vpos[vpos_p] = vp;
			vval[vpos_p] = drand48();
			vp++; vpos_p++;
		}
	}

	fp_t *x = malloc(DIM*sizeof(fp_t));
	for(int i = 0; i < DIM; ++i)
		x[i] = drand48();
	fp_t *y = malloc(DIM*sizeof(fp_t));
#ifdef SINGLE_PRECISION
	mkl_scsrmv("N", &m, &m, &alpha, "GLNC", vval, vpos, vptr, vptr+1, x, &beta, y);
	mkl_scsrsv("N", &m, &alpha, "TLNC", vval, vpos, vptr, vptr+1, x, y);
	mkl_cspblas_scsrgemv("N", &m, vval, vptr, vpos, x, y);
#else
	mkl_dcsrmv("N", &m, &m, &alpha, "GLNC", vval, vpos, vptr, vptr+1, x, &beta, y);
	mkl_dcsrsv("N", &m, &alpha, "TLNC", vval, vpos, vptr, vptr+1, x, y);
	mkl_cspblas_dcsrgemv("N", &m, vval, vptr, vpos, x, y);
#endif
	
	free(x); free(y);
	free(vptr); free(vpos); free(vval);
	free(t);
	printf("warm-up passed\n");
#endif
}
