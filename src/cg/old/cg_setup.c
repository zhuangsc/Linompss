#include "cg_setup.h"

#include "fptype.h"
#include "fpblas.h"
#include "async_struct.h"
#include "genmat.h"
#include "convertmat.h"
#include "matfprint.h"
#include "matfread.h"
#include "ompss_apps.h"
#include "genmat_config.h"
#include "async_struct.h"

#include <stdlib.h>
#include <stdio.h>

#if USE_SPARSE
#include "hb.h"
#include "hbconvrt.h"
#endif


extern unsigned long works;
extern char *rhsfname;
extern void *Ahhb;
static fp_t lmarker[4];


typedef struct {
	async_t sync[3]; 
	fp_t bstart;
} worklay_t;


int cg_setup(int n, int bm, int bn, int s, void **A, fp_t **x, fp_t **rhs, fp_t **xstar, fp_t **work) 
{
#if USE_SPARSE
	hbmat_t *lA = *A;
	Ahhb = hb2hbh(lA, bm, 1);
//	*A = hb2hbh_symcsr(lA, bm, NULL, 0);
#else
	fp_t *lA = *A = GENMAT_CONFIG(ompssapp_CG, "A", n, n, n, GENMATCONF_COVAR);
	if ( lA == NULL ) {
		fprintf(stderr, "setup: could not allocate A\n");
		return 5;
	}
#endif

	int ns = n * s;
	int dupl = 4;

	fp_t *lx = *x = malloc(dupl * ns * sizeof(fp_t));
	if ( lx == NULL ) {
		fprintf(stderr, "err: could not allocate x\n");
		return 1;
	}
	

/* Construct or read the rhs (B) */
	fp_t *lrhs;
	if ( rhsfname != NULL ) {
		lrhs = *rhs = (fp_t*) calloc(ns, sizeof(fp_t));
		if ( lrhs == NULL ) {
			fprintf(stderr, "err: could not allocate rhs\n");
			return 3;
		}

		FILE *rhsf;
		if ( rhsf = fopen(rhsfname, "r") ) {
			READ_MM2DENSE(rhsf, n, s, lrhs);
			fclose(rhsf);
		} else {
			fprintf(stderr, "setup: could not open %s\n", rhsfname);
			return 6;
		}
	} else {
#if 0
		fp_t *lxstar = *xstar = GENMAT(n, s, n);
		if ( lxstar == NULL ) {
			fprintf(stderr, "err: could not allocate xstar\n");
			return 2;
		}
		FPRINT_DENSE2MM("xstar.mm", "xstar", n, s, lxstar, n);

		BLAS_gemm(OMPSSBLAS_TRANSP, OMPSSBLAS_NTRANSP, n, s, n, FP_ONE, lA, n, lxstar, n, FP_NOUGHT, lrhs, n); 
		// or rather mkl_dcsrsymv ("U", &n, mat.vval, mat.vptr, mat.vpos, x, rhs);
#endif
		lrhs = *rhs = GENMAT_CONFIG(ompssapp_CG, "rhs", n, s, n, GENMATCONF_COSTAS);
		*xstar = NULL;
	}

/* align work buffer on 128-byte boundary */
	int nsalgn = ((ns + 0x7f ) >> 7 ) << 7;
	int salgn = ((s + 0x7f ) >> 7 ) << 7;
	int tmps = nsalgn * dupl * 4 + 3 * salgn * dupl;

	works = sizeof(worklay_t) + (tmps - 1) * sizeof(fp_t);  
	/* an additional 4 for the sentinels */
	if ( posix_memalign(work, 128, works + 4 * sizeof(fp_t)) ) {
		fprintf(stderr, "setup: failed to allocate aligned memory\n");
		return 4;
	}

	worklay_t *worklay = (worklay_t*) *work;
	async_setup(worklay->sync);

	fp_t *marker = &worklay->bstart + tmps;
	marker[0] = lmarker[0] = FP_RAND();
	marker[1] = lmarker[1] = FP_RAND();
	marker[2] = lmarker[2] = FP_RAND();
	marker[3] = lmarker[3] = FP_RAND();

	works = tmps;

	return 0;
}


int cg_cleanup(int n, int s, void *A, fp_t *x, fp_t *xstar, fp_t *rhs, fp_t *work) 
{
	worklay_t *worklay = (worklay_t*) work;

	async_t *sync = worklay->sync;
	async_fini(sync, 3);

	fp_t *marker = &worklay->bstart + works;
	int ok = marker[0] == lmarker[0] && marker[1] == lmarker[1] && marker[2] == lmarker[2] && marker[3] == lmarker[3];
	if ( ! ok ) {
		fprintf(stderr, "err: marker does not look good\n");
	}

	free(A);
	free(x);
	free(xstar);
	free(rhs);
	free(work);
}
