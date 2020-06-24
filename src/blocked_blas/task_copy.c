#include <omp.h>
#include "task_copy.h"

#include "fptype.h"
#include "fpblas.h"
#include "blas.h"
#include "selfsched.h"


#ifdef DOUBLE_PRECISION

#define __t_copy 					task_dcopy
#define __t_copy_sched 				task_dcopy_sched

#else

#define __t_copy 					task_scopy
#define __t_copy_sched 				task_scopy_sched

#endif


void __t_copy(int p, int bm, int bn, int m, int n, fp_t *X, fp_t *Y) 
{
	int j;
	for ( j=0; j<bn; ++j ) {
		BLAS_cp(bm, &X[j*m], i_one, &Y[j*m], i_one);
	}
}


#if 0
void __t_copy_sched(async_t *sync, int z, int p, int bm, int bn, int m, int n, fp_t *X, int iy, selfsched_t *schedY, fp_t *Y) 
{
	int j;
	for ( j=0; j<bn; ++j ) {
		BLAS_cp(bm, &X[j*m], i_one, &Y[j*m], i_one);
	}

	fp_t len = BLAS_dot(bm, X, I_ONE, X, I_ONE);

	pthread_mutex_t *mutex = &sync->mutex;
	pthread_mutex_lock(mutex);

	//prof_add(&sync->prof, iy, len);
	//sched_done(schedY, z, iy);

	pthread_mutex_unlock(mutex);
}
#endif
