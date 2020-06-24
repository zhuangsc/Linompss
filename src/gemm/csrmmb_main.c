#include "csrmmb_main.h"

//#include "bblas_gemm.h"
#include "fptype.h"
#include "hb.h"
#include "bsblas_csrmmb.h"


/* d and c not supported yet */
int CSRMMB(int b, int d, int c, int n, fp_t alpha, hbmat_t *Ah, fp_t *B, int ldb, fp_t beta, fp_t *C, int ldc) 
{

	BSBLAS_CSRMMB(b, n, alpha, Ah, B, ldb, beta, C, ldc); 

	return 0;
}
