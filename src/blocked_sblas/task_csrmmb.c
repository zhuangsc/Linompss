#include "task_csrmmb.h"


#include "fptype.h"
#include "fpsblas.h"
#include "matfprint.h"
#include "selfsched.h"
#include "hbconvrt.h"


#ifdef SINGLE_PRECISION

#define _t_csrmmb		task_scsrmmb

#else

#define _t_csrmmb		task_dcsrmmb

#endif


void _t_csrmmb(int n, fp_t alpha, char *descra, hbmat_t *A, fp_t *B, int ldb, fp_t beta, fp_t *C, int ldc, int si, selfsched_t *sched) 
{
	int m = A->m;
	int k = A->n;

	/* block has not yet been created */
	if ( A->vval == NULL ) {
		hbmat_t *Ahb = A->orig;
		hbmat_t *Ahbh = A->hyper;
		int I = A->orig_row;
		int J = A->orig_col;
		int b = Ahbh->b;

		HB2HBH_block(I, J, Ahb, b, A);	
	}

	int *vpos = A->vpos;
	int *vptr = A->vptr;
	int *vval = A->vval;

	if ( sched != NULL ) {
		sched_wait(sched, si);
	}

#if 0
	printf("\n\n\n");
	PRINT_HB(stdout, "A", A, 1);
	PRINT_DENSE2MM(stdout, "B", k, n, B, ldb);
	PRINT_DENSE2MM(stdout, "C", m, n, C, ldc);

	printf("%i %i %i | B %p %i C %p %i\n", m, n, k, B, ldb, C, ldc);
#endif
	SBLAS_csrmm("N", m, n, k, alpha, descra, vval, vpos, vptr, vptr+1, B, ldb, beta, C, ldc);
	//PRINT_DENSE2MM(stdout, "C", m, n, C, ldc);
	//fflush(0);
}

