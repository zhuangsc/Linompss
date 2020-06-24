#include "ompss_sparse_gs.h"

#include <stdio.h>
#include "hb.h"
#include "hbext.h"
#include "hbconvrt.h"
#include "gauss-seidel_main.h"



/* sparse, 0-based hbmat_t */
void ompss_csr_gs(int bs, hbmat_t *A,  fp_t *v_x, fp_t *v_b, int max_iter, int lookahead, fp_t threshold, fp_t *work) 
{

	int dim = A->n;
	hbmat_t *Ahbh = hb2hbh(A, bs, 1);

	int *vdiag = Ahbh->vdiag;
	int N = Ahbh->n;
	hbmat_t **diagL = malloc ( N * sizeof(hbmat_t*));
	hbmat_t **A1 = malloc ( N * sizeof(hbmat_t*));
	hbmat_t **vval = Ahbh->vval;
	int I;
	for ( I = 0; I < N; ++I ){
		int *etree0 = etree(vval[vdiag[I]]);
		A1[I] = hb2hbh_symcsr(vval[vdiag[I]], vval[vdiag[I]]->m, etree0, 1);
		diagL[I] = ((hbmat_t**)A1[I]->vval)[0];
		diagL[I]->FACT = 0;
	}
	GS_MAIN_CSR(Ahbh, v_x, v_b, bs, max_iter, diagL, lookahead, threshold, work);

}
