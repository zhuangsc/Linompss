#include "ompss_sparse_chol.h"

#include "chols_llmain.h"

#include <stdio.h>

/* sparse, 0-based hbmat_t */
hbmat_t* ompss_csr_chol_ll(int b, hbmat_t *A, int *work) 
{

	int *tree = etree(A);
	hbmat_t *Acsr = hb2hbh_hyper_sym_csr(A, b, tree);

	hbmat_t *Acsc = Acsr->trans = csc2csr(Acsr);

	chols_ll_upper(Acsr, Acsc, work);

	return Acsr;

}

hbmat_t* ompss_csc_chol_ll(int b, hbmat_t *A, int *work) 
{

	hbmat_t *tAcsr = A->trans = csc2csr(A);
	int *tree = etree(tAcsr);
	free(tAcsr);

	hbmat_t *Acsc = hb2hbh_sym_etree(A, b, tree);
	hbmat_t *Acsr = Acsc->trans = csc2csr(Acsc);

	free(tree);

	chols_ll(Acsc, Acsr, work);

	return Acsc;

}
