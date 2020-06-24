#include <time.h>
#include <sys/time.h>

#include "jacobi_check.h"
#include "bsblas_gemv_csr.h"
#include "densutil.h"

extern fp_t *v_b;
extern fp_t *v_x0;
extern int dim;
extern hbmat_t *Ahbh;
extern hbmat_t *Ahb;
extern int bs;

fp_t jacobi_check(int check, int res_p)
{
	int res = 0;
	if ( check ) {
		fp_t x0 = 0;
		fp_t *v_res = &v_x0[res_p];
//		fp_t *v_res = iter % 2 ? (&v_x0[dim]) : v_x0; //The position the latest result resides depends on the number of iterations
		fp_t *cnorm_x0 = malloc( dim * sizeof(fp_t));
		int i;
		for ( i = 0; i < dim; ++i )
			cnorm_x0[i] = v_b[i];
		fp_t alpha = FP_MONE; fp_t beta = FP_ONE;
		SBLAS_csrmv("N", &(Ahb->m), &(Ahb->n), &alpha, "GLNC", Ahb->vval, Ahb->vpos, Ahb->vptr, (Ahb->vptr)+1, v_res, &beta, cnorm_x0);
		x0 = VECTOR_2NORM(cnorm_x0, dim);
		fp_t norm_b = VECTOR_2NORM(v_b, dim);
		printf("2-norm B-A*x0: %e\n", x0/norm_b);
		res = x0/norm_b;
	}
	return res;
}
