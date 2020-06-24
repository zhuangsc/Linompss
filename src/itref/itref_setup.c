#include "itref_setup.h"

#include "fptype.h"
#include "fpblas.h"
#include "genmat.h"
#include "async_struct.h"
#include "convertmat.h"
#include "matfprint.h"
#include "matfread.h"
#include "ompss_apps.h"
#include "genmat_config.h"

#include <stdlib.h>
#include <stdio.h>

#if USE_SPARSE
#include "hb.h"
#include "hbconvrt.h"
#include "dcsparse.h"
extern void *Ahb;
extern void *Acsr;
extern double *dvval;
extern void *preconditioner;
extern css **S;
extern csn **N;
extern cs *Acs;
#endif

extern unsigned long works;
extern char *rhsfname;
extern char *aname;
static double lmarker[4];
extern double *mat_energy;


typedef struct {
	async_t sync[3]; 
	fp_t bstart;
} worklay_t;

void dgen_cg_rhs(int m, int n, double *b)
{
	unsigned int seed = 1591613054;
	srand(seed);
	double range = (double) 1;

	int i;
	for ( i = 0; i < m*n; i++ ) {
		b[i] = ((double)rand() / (double)RAND_MAX ) * range - range/2;
	}
}

int itref_dsetup(int n, int bm, int bn, int s, void **A, double *x[], double **rhs, double **work) 
{
#if USE_SPARSE
	hbmat_t *lAhb = Acsr;
	int elemc = lAhb->elemc;
	dvval = lAhb->vval;
//	hbmat_t *lAhbs = Acsrs = malloc(sizeof(hbmat_t));
//	hb_init_basic(Acsrs, Acsr);
//	lAhbs->vptr = lAhb->vptr; lAhbs->vpos = lAhb->vpos;
//	svval = malloc(elemc * sizeof(double));
//	lAhbs->vval = svval;
//	int i;
//	for ( i = 0; i < elemc; i++ ) {
//		svval[i] = dvval[i];
//	}
	mat_energy = hb_energy_row_double(Acsr, bm);
	int bs = (n+bm-1)/bm;
	preconditioner = malloc(bs * sizeof(hbmat_t));
	hb_sym_diag_block_double(Acsr, bm, preconditioner);
	hbmat_t *precond = preconditioner;
	Acs = malloc(bs*sizeof(cs));
	S = malloc(bs*sizeof(css*));
	N = malloc(bs*sizeof(csn*));
	for (int i = 0; i < bs; i++ ) {
		hbmat_t *Apre = &precond[i];
		cs *Acspre = &Acs[i];
		Acspre->nz = -1;
		Acspre->nzmax = Apre->elemc;
		Acspre->n = Acspre->m = Apre->m;
		Acspre->p = Apre->vptr;
		Acspre->i = Apre->vpos;
		Acspre->x = Apre->vval;
		S[i] = cs_schol_double(Acspre, 0);
		N[i] = cs_chol_double(Acspre, S[i]);
	}
	for (int i = 0; i < bs; i++) {
		free(precond[i].vptr);
		free(precond[i].vpos);
		free(precond[i].vval);
	}
	free(precond);

#else
	double *lA;
	double *fA;
	if ( aname != NULL ) {
		printf("read A\n");
		lA = *A = (double*) calloc(n*n, sizeof(double));
		fA = calloc(n*n, sizeof(double));
		FILE *laf;
		if ( laf = fopen(aname, "r")) {
			READ_MM2DENSE(laf, n*n, 1, fA);
			fclose(laf);
			int i;
			for (i=0; i<n*n; ++i)
				lA[i] = (double) fA[i];
		} else {
			fprintf(stderr, "setup: could not open %s\n", aname);
			return 7;
		}
	} else {
		printf("gen A\n");
		lA = *A = dgenmat_config(ompssapp_ITREF, "A", n, n, n, GENMATCONF_COVAR);
		int i;
		fA = calloc(n*n, sizeof(double));
		for (i=0; i<n*n; ++i)
			fA[i] = (double) lA[i];
		FPRINT_DENSE2MM("A.dat", "A", n, n, fA, n);
	}
	if ( lA == NULL ) {
		fprintf(stderr, "setup: could not allocate A\n");
		return 5;
	}

	//TODO Matrix energy distribution
	int szeb = (n+bm-1)/bm;
	mat_energy = calloc(szeb+1, sizeof(double));
	mat_energy[szeb] = dmat_energy_sym(n, n, bm, lA, mat_energy);
//	double *tmpe = calloc(szeb+1, sizeof(double));
//	tmpe[szeb] = dmat_energy_sym(n, n, bm, lA, tmpe);
//	int i;
//	for ( i = 0; i < szeb+1; i++ ) {
//		mat_energy[i] = tmpe[i];
//	}
//	free(tmpe);

#endif

	int ns = n * s;
	int dupl = 4;

	double *lx = malloc(2 * dupl * ns * sizeof(double));
	if ( lx == NULL ) {
		fprintf(stderr, "err: could not allocate x\n");
		return 1;
	}
	x[0] = &lx[0];
	x[1] = &lx[dupl*ns];

	/* Construct or read the rhs (B) */
	double *lrhs;
	double *frhs;
	if ( rhsfname != NULL ) {
		printf("read rhs\n");
		lrhs = *rhs = (double*) calloc(ns, sizeof(double));
		frhs = calloc(ns, sizeof(double));
		if ( lrhs == NULL ) {
			fprintf(stderr, "err: could not allocate rhs\n");
			return 3;
		}

		FILE *rhsf;
		if ( rhsf = fopen(rhsfname, "r") ) {
			READ_MM2DENSE(rhsf, n*s, 1, frhs);
			fclose(rhsf);
			int i;
			for(i=0; i< n*s; ++i)
				lrhs[i] = (double) frhs[i];
		} else {
			fprintf(stderr, "setup: could not open %s\n", rhsfname);
			return 6;
		}
	} else {
		printf("gen rhs\n");
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
#if USE_SPARSE
		lrhs = *rhs = malloc(n*s*sizeof(double));
		dgen_cg_rhs(n,s,lrhs);
#else
		lrhs = *rhs = dgenmat_config(ompssapp_ITREF, "rhs", n, s, n, GENMATCONF_RAND);
#endif
		frhs = calloc(n*s, sizeof(double));
		int i;
		for ( i = 0; i < n*s; ++i )
			frhs[i] = lrhs[i];
		FPRINT_DENSE2MM("RHS.dat", "RHS", n, s, frhs, n);
	}

	/* align work buffer on 128-byte boundary */
	int nsalgn = ((ns + 0x7f ) >> 7 ) << 7;
	int salgn = ((s + 0x7f ) >> 7 ) << 7;
	int tmps = nsalgn * dupl * 4 + 3 * salgn * dupl;

	works = sizeof(worklay_t) + (tmps - 1) * sizeof(fp_t);  
	/* an additional 4 for the sentinels */
	if ( posix_memalign(work, 128, works + 4 * sizeof(double)) ) {
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


int itref_cleanup(int n, int s, void *A, double *x[], double *rhs, float *work)
{
	worklay_t *worklay = (worklay_t*) work;

	async_t *sync = worklay->sync;
	async_fini(sync, 3);

	fp_t *marker = &worklay->bstart + works;
	int ok = marker[0] == lmarker[0] && marker[1] == lmarker[1] && marker[2] == lmarker[2] && marker[3] == lmarker[3];
	if ( ! ok ) {
		fprintf(stderr, "err: marker does not look good\n");
	}

#ifdef USE_SPARSE
	hb_free(A);
#else
	free(A);
#endif
	free(x[0]); 
	free(rhs);
	free(work);
}
