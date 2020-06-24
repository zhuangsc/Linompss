#include "hb.h"


void hb_free(hbmat_t *A)
{
//	if ( A->b != 0 ) {
//		int elemc = A->elemc;
//		hbmat_t **vval = A->vval;
//		int c;
//		for( c=0; c<elemc; ++c ) {
//			hb_free(vval[c]);
//		}
//	}

	free(A->vptr); 
	free(A->vpos);
	free(A->vval);
	free(A->vdiag);
	free(A->e_tree);

	free(A);
}

void hbh_free2(hbmat_t *A)
{
	pthread_mutex_destroy(A->mtx);
	int elemc = A->elemc;
	free(A->e_tree);
	free(A->vptr);
	free(A->vpos);
	free(((hbmat_t**)A->vval)[0]);
	free(A->vptr_pool);
	free(A->vpos_pool);
	free(A->vval_pool);
//	hb_free(A);

}
