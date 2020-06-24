#include "cg_setup.h"


#include "fptype.h"
#include "fpblas.h"
#include "blas.h"
#include "hb.h"
#include "hbconvrt.h"
#include "genmat.h"
#include "matfprint.h"

#include <stdlib.h>
#include <stdio.h>


extern unsigned long works;
extern char *rhsfname;
static fp_t lmarker[4];
static fp_t **lwork;

int cg_setup(int n, int bm, int bn, int s, void **A, fp_t **x, fp_t **rhs, fp_t **xstar, fp_t **work) {
	hbmat_t *lA = *A;;
	//one2zero(Ahb);
	*A = hb2hbh_symcsr(lA, bm, NULL, 0);

	//PRINT_HB(stdout, "A", *A, 0);

	int ns = n * s;
	int dupl = 4;

	fp_t *lx = *x = malloc(dupl * ns * sizeof(fp_t));
	if ( lx == NULL ) {
		fprintf(stderr, "err: could not allocate x\n");
		return 1;
	}
	

/* Construct or read the rhs (B) */
#if 0
	fp_t *lrhs = *rhs = (fp_t*) calloc(ns, sizeof(fp_t));
	if ( lrhs == NULL ) {
		fprintf(stderr, "err: could not allocate rhs\n");
		return 3;
	}

	if ( rhsfname != NULL ) {
		*xstar = NULL;
		FILE *rhsf;
		if ( rhsf = fopen(rhsfname, "r") ) {
			READ_MM2DENSE(rhsf, n, s, lrhs);
			fclose(rhsf);
		} else {
			fprintf(stderr, "setup: could not open %s\n", rhsfname);
			return 6;
		}
	} else {
		fp_t *lxstar = *xstar = GENMAT(n, s);
		if ( lxstar == NULL ) {
			fprintf(stderr, "err: could not allocate xstar\n");
			return 2;
		}
		FPRINT_DENSE2MM("xstar.mm", "xstar", n, s, lxstar);

		BLAS_gemm(MAT_NOTRANSP, MAT_NOTRANSP, n, s, n, FP_ONE, lA, n, lxstar, n, FP_NOUGHT, lrhs, n); 
		FPRINT_DENSE2MM("rhs.mm", "rhs", n, s, lrhs);
	}
#else
		fp_t *lrhs = *rhs = GENMAT_COSTAS(n, s);
		FPRINT_DENSE2MM("rhs.mm", "rhs", n, s, lrhs, n);
		*xstar = NULL;
#endif

/* align work buffer on 128-byte boundary */
	int nsalgn = ((ns + 0x7f ) >> 7 ) << 7;
	int salgn = ((s + 0x7f ) >> 7 ) << 7;
	works = nsalgn * (dupl * 3) + 2 * salgn * dupl;
	/* an additional 4 for the sentinels */
	if ( posix_memalign(work, 128, (works + 4 )* sizeof(fp_t)) ) {
		fprintf(stderr, "setup: failed to allocate aligned memory\n");
		return 4;
	}

	fp_t *marker = *work + works;
	marker[0] = lmarker[0] = FP_RAND();
	marker[1] = lmarker[1] = FP_RAND();
	marker[2] = lmarker[2] = FP_RAND();
	marker[3] = lmarker[3] = FP_RAND();

	return 0;
}

int cg_setup_check() {
}


int cg_cleanup(int n, int s, void *A, fp_t *x, fp_t *xstar, fp_t *rhs, fp_t *work) {
	fp_t *marker = work + works;
	int ok = marker[0] == lmarker[0] && marker[1] == lmarker[1] && marker[2] == lmarker[2] && marker[3] == lmarker[3];
	if ( ! ok ) {
		fprintf(stderr, "err: marker does not look good\n");
	}

	hb_free(A);
	free(x);
	free(xstar);
	free(rhs);
	free(work);
}
