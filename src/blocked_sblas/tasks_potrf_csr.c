#include "tasks_potrf_csr.h"


#include "symfac.h"
#include "blas.h"
#include "fpblas.h"


#ifdef DOUBLE_PRECISION

#define __t_potrf_csr		task_dpotrf_csr

#else

#define __t_potrf_csr		task_spotrf_csr

#endif


void __t_potrf_csr(hbmat_t* A) 
{
 	/* Check if the input matrix has properly set */
	if ( A->vptr == NULL ) {
		hyper_sym_csr_task2(A);
	}

	int m = A->m;
	int* vptr = A->vptr; int* vpos = A->vpos; 
	fp_t* vval = A->vval;

	int d;
	for ( d=0 ; d<m; ++d ) {
		int dstart = vptr[d];
		int dend = vptr[d+1];

		int diagidx = dend - 1;
		fp_t diagel = vval[diagidx];
		/* compute element on diagonal (d,d) */
		int j;
		for ( j = dstart; j < diagidx; ++j ) {
			diagel -= vval[j]*vval[j];
		}
		diagel = sqrt(diagel);
		vval[diagidx] = diagel;


		int i;
		for ( i = d+1; i < m; ++i ) {
			int jx = vptr[i];
			int j = vpos[jx];
			int jend = vptr[i+1];

			int djx = dstart;
			int dj = vpos[djx];

			fp_t acc = 0;
			while ( dj < d && djx < dend && j < d && jx < jend ) {
				int djadv = dj < j; 
				int jadv  = j  < dj; 

				jx  += jadv? 1 : 0;
				djx += djadv? 1 : 0;
				dj = vpos[djx];
			 	j  = vpos[ jx];

				if ( dj == j && j < d && djx < dend && jx < jend ) {
					acc += vval[djx++] * vval[jx++];
					dj = vpos[djx];
			 		j  = vpos[ jx];
				}
			}

			while ( vpos[jx] < d && jx < jend ) { ++jx; }

			if ( jx < jend && vpos[jx]==d ) {
				fp_t subdval = vval[jx] - acc;
				vval[jx] = subdval / diagel;
			}
		}
	}
}
