#include "cg_setup.h"

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

#include "hb.h"
#include "hbconvrt.h"
#include "dcsparse.h"

extern void *Ahbh;
extern void *Ahb;
extern void *Acsr;
extern void *preconditioner;
extern css **S;
extern csn **N;
extern cs *Acs;
extern char *rhsfname;
extern char *aname;


void gen_cg_rhs(int m, int n, double *b)
{
	unsigned int seed = 1591613054;
	srand(seed);
	double range = (double) 1;

	int i;
	for ( i = 0; i < m*n; i++ ) {
		b[i] = ((double)rand() / (double)RAND_MAX ) * range - range/2;
	}
}

int cg_setup(int n, int bm, double **x, double **rhs) 
{
	Ahbh = hb2hbh_double(Acsr, bm, 1);
	hbmat_t *lAhb = Acsr;
	int elemc = lAhb->elemc;
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

	*x = malloc(n * sizeof(double));
//	double *lx = malloc(2 * n * sizeof(double));
//	if ( lx == NULL ) {
//		fprintf(stderr, "err: could not allocate x\n");
//		return 1;
//	}
//	x[0] = &lx[0];
//	x[1] = &lx[n];

	/* Construct or read the rhs (B) */
	
	double *lrhs;
	if ( rhsfname != NULL ) {
		printf("read rhs\n");
		lrhs = *rhs = (double*) calloc(n, sizeof(double));
		if ( lrhs == NULL ) {
			fprintf(stderr, "err: could not allocate rhs\n");
			return 3;
		}

		FILE *rhsf;
		if ( rhsf = fopen(rhsfname, "r") ) {
			READ_MM2DENSE(rhsf, n, 1, lrhs);
			fclose(rhsf);
		} else {
			fprintf(stderr, "setup: could not open %s\n", rhsfname);
			return 6;
		}
	} else {
		printf("gen rhs\n");

		lrhs = *rhs = malloc(n*sizeof(double));
		gen_cg_rhs(n, 1, lrhs);
		FPRINT_DENSE2MM("RHS.dat", "RHS", n, 1, lrhs, n);
	}

	return 0;

}
