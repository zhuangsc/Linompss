#include "csrmm_main.h"

#include "fptype.h"
#include "hb.h"
#include "bsblas_csrmm.h"
#include "selfsched.h"


/* d and c not supported yet */
int CSRMM(int b, int d, int c, int n, fp_t alpha, hbmat_t *Ahbh, fp_t *B, int ldb, fp_t beta, fp_t *C, int ldc, selfsched_t *sched) 
{
	/* n is the number of of cols of B or C */
	BSBLAS_CSRMM(b, n, alpha, Ahbh, B, ldb, beta, C, ldc, sched); 

	return 0;
}
