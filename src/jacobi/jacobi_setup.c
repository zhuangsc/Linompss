#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include "jacobi_setup.h"

#include "fpsblas.h"
#include "fptype.h"
#include "hb.h"
#include "hbext.h"
#include "hbconvrt.h"
#include "iohb.h"
#include "genmat.h"
#include "symfac.h"



extern void *Ahb;
extern void *Ahbh;
extern fp_t *v_b;
extern fp_t *v_x;
extern fp_t *v_x0;
extern int format;
extern int bs;
extern int dim;
extern hbmat_t **diagL;
extern hbmat_t **A1;
extern fp_t *work;

int jacobi_setup (const char *fname) 
{
	Ahb = (hbmat_t*) malloc(sizeof(hbmat_t));
	hbmat_t *A = Ahb;
	hb_reset(Ahb);

	readHB_newmat(fname, &(A->m), &(A->n), &(A->elemc), &(A->vptr), &(A->vpos), (fp_t **)&(A->vval));

	one2zero(A);

	dim = A->n ;
	Ahbh = hb2hbh(A, bs, format);
	hbmat_t *A0 = Ahbh;
	v_b = malloc(dim * sizeof(fp_t));
	v_x0 = calloc(2 * dim, sizeof(fp_t));
	work = malloc((1+2*dim) * sizeof(fp_t));

	v_x = GENMAT(dim, 1, dim);
	fp_t alpha = FP_ONE; fp_t beta = FP_NOUGHT;
	SBLAS_csrmv("N", &A->m, &A->n, &alpha, "GLNC", A->vval, A->vpos, A->vptr, (A->vptr)+1, v_x, &beta, v_b);

	int *vdiag = A0->vdiag;
	int N = A0->n;
	diagL = malloc ( N * sizeof(hbmat_t*));
	A1 = malloc ( N * sizeof(hbmat_t*));
	hbmat_t **vval = A0->vval;
	int I;
	for ( I = 0; I < N; ++I ){
		int *etree0 = etree(vval[vdiag[I]]);
		A1[I] = hb2hbh_symcsr(vval[vdiag[I]], vval[vdiag[I]]->m, etree0, 1);
		diagL[I] = ((hbmat_t**)A1[I]->vval)[0];
		diagL[I]->FACT = 0;
	}
	return 0;
}

void jacobi_shutdown () 
{
	hb_free(Ahb);
	hb_free(Ahbh);
	int N = dim/bs;
	int I;
	for ( I = 0; I < N; ++I )
		free(A1[I]);
	free(v_b); free(v_x); free(v_x0);
	free(work);
}
