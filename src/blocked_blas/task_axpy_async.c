#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "task_axpy_async.h"

#include "fptype.h"
#include "fpblas.h"
#include "blas.h"
#include "async_struct.h"
#include "task_log.h"




#ifdef DOUBLE_PRECISION

#define __t_axpy_async				task_daxpy_async
#define __t_cpaxpy_async			task_dcpaxpy_async
#define __t_scal_axpy_async			task_scal_daxpy_async
#define __t_scal_cpaxpy_async		task_scal_dcpaxpy_async
#define __t_scal_cpaxpy_comb_async	task_scal_dcpaxpy_comb_async
#define __t_scal_cpaxpy_comb_async_release task_scal_dcpaxpy_comb_async_release
#define __t_scal_cpaxpy_comb_async_concurrent	task_scal_dcpaxpy_comb_async_concurrent
#define __t_cpaxpy_comb4_async		task_dcpaxpy_comb4_async

#else

#define __t_axpy_async				task_saxpy_async
#define __t_cpaxpy_async			task_scpaxpy_async
#define __t_scal_axpy_async			task_scal_saxpy_async
#define __t_scal_cpaxpy_async		task_scal_scpaxpy_async
#define __t_scal_cpaxpy_comb_async	task_scal_scpaxpy_comb_async
#define __t_scal_cpaxpy_comb_async_release task_scal_scpaxpy_comb_async_release
#define __t_scal_cpaxpy_comb_async_concurrent	task_scal_scpaxpy_comb_async_concurrent
#define __t_cpaxpy_comb4_async		task_scpaxpy_comb4_async

#endif


/* unused */
void __t_axpy_async(int dotid, async_t *sync, int bm, int bn, int m, int n, fp_t *Anum, fp_t *Aden, fp_t *X, fp_t *D, fp_t *Y) {
	pthread_mutex_t *mutex = &sync->mutex;
	pthread_mutex_lock(mutex);
	if ( sync->create != dotid ) {
		++sync->wait;		
		pthread_cond_wait(&sync->cond, mutex);
	} else if ( sync->consume != dotid ) {
		sync->consume = dotid;

	}
	pthread_mutex_unlock(mutex);

	int j;
	for ( j=0; j<bn; ++j) {
		fp_t factor = Anum[j] / Aden[j];
		BLAS_axpy(bm, factor, X, i_one, Y, i_one);
		X += m;
		Y += m;
	}
}

void __t_scal_axpy_async(int dotid, async_t *sync, int bm, int bn, int m, int n, fp_t alpha, fp_t *Anum, fp_t *Aden, fp_t *X, fp_t *Y) {
	int log = 0;	
	pthread_mutex_t *mutex = &sync->mutex;
	int pcnt;

	pthread_mutex_lock(mutex);
	if ( sync->create != dotid ) {
		++sync->wait;		
		pthread_cond_wait(&sync->cond, mutex);
	} else if ( sync->consume != dotid ) {
		sync->consume = dotid;
		log = 1;
	}
	pcnt = sync->pcnt;
	pthread_mutex_unlock(mutex);

#if 0
	// scale if enabled and reset
	int scale = sync->flags & ASYNC_SCALE;
	fp_t corr = ((fp_t)sync->pcompl) / ((fp_t)pcnt);
	corr = scale ? corr : FP_ONE;
#endif


	fp_t localA[bn];
	int j;
	for ( j = 0; j < bn; ++j ) {	
		localA[j] = alpha * Anum[j] / Aden[j];
	}

	for ( j=0; j<bn; ++j) {
		BLAS_axpy(bm, localA[j], X, i_one, Y, i_one);
		X += m;
		Y += m;
	}

#if 0
	if ( log && sync->log ) {
		async_log(sync, pcnt, corr);
	}
#endif
}


/* Z = Z + Anum / Aden * X */
void __t_cpaxpy_async(int dotid, async_t *sync, int bm, int bn, int m, int n, fp_t *Anum, fp_t *Aden, fp_t *X, fp_t *D, fp_t *Y, fp_t *Z) {
	int log = 0;
	pthread_mutex_t *mutex = &sync->mutex;
	int pcnt;

	pthread_mutex_lock(mutex);
	if ( sync->create != dotid ) {
		++sync->wait;		
		pthread_cond_wait(&sync->cond, mutex);
	} else if ( sync->consume != dotid ) {
		sync->consume = dotid;
		log = 1;
	}
	pcnt = sync->pcnt;
	pthread_mutex_unlock(mutex);

#if 0
	// scale if enabled and reset
	int scale = sync->flags & ASYNC_SCALE;
	fp_t corr = ((fp_t)sync->pcompl) / ((fp_t)pcnt);
	corr = scale ? corr : FP_ONE;
#endif


	int j;
	for ( j=0; j<bn; ++j) {
		fp_t factor = Anum[j] / (Aden[j]  * 1.0);
		BLAS_cp(bm, Y, i_one, Z, i_one);
		BLAS_axpy(bm, factor, X, i_one, Z, i_one);
		X += m;
		Y += m;
		Z += m;
	}

#if 0
	if ( log && sync->log ) {
		async_log(sync, pcnt, corr);
	}
#endif
}


void __t_scal_cpaxpy_async(int dotid, async_t *sync, int bm, int bn, int m, int n, fp_t alpha, fp_t *Anum, fp_t *Aden, fp_t *X, fp_t *Y, fp_t *Z) {
	int pcnt; 
	pthread_mutex_t *mutex = &sync->mutex;
	int log = 0;

	pthread_mutex_lock(mutex);
	if ( sync->create != dotid ) {
		++sync->wait;		
		pthread_cond_wait(&sync->cond, mutex);
	} else if ( sync->consume != dotid ) {
		sync->consume = dotid;
		log = 1;
	}
	pcnt = sync->pcnt;
	pthread_mutex_unlock(mutex);

#if 0
	// scale if enabled and reset
	int scale = sync->flags & ASYNC_SCALE;
	fp_t corr = ((fp_t)sync->pcompl) / ((fp_t)pcnt);
	corr = scale ? corr : FP_ONE;
#endif

	fp_t localA[bn];
	int j;
	for ( j = 0; j < bn; ++j ) {	
		localA[j] = alpha * Anum[j] / (Aden[j] * 1.0);
	}

	for ( j=0; j<bn; ++j) {
		BLAS_cp(bm, Y, i_one, Z, i_one);
		BLAS_axpy(bm, localA[j], X, i_one, Z, i_one);
		X += m;
		Y += m;
	}

#if 0
	/* logging */
	if ( log && sync->log ) {
		async_log(sync, pcnt, corr);
	}
#endif
}

void __t_scal_cpaxpy_comb_async(int dotid, int l, async_t *sync, int p, int bm, int bn, int m, int n, 
	fp_t alpha, fp_t *Anum, fp_t *Aden, fp_t *X1, fp_t *X2, fp_t *Y1, fp_t *Y2, fp_t *Z1, fp_t *Z2) 
{
	int pcnt; 
	pthread_mutex_t *mutex = &sync->mutex;
	int log = 0;

	BLAS_cp(bm, Y2, i_one, Z2, i_one);
	fp_t len = BLAS_dot(bm, X2, I_ONE, Y1, I_ONE);

	pthread_mutex_lock(mutex);
	if ( sync->create != dotid ) {
		++sync->wait;		
		pthread_cond_wait(&sync->cond, mutex);
	} else if ( sync->consume != dotid ) {
		sync->consume = dotid;

		pcnt = sync->pcnt;
		log = 1;
	}
	pthread_mutex_unlock(mutex);


#if 1
	int j;
	for ( j=0; j<bn; ++j) {
		/* update of x */
		fp_t factor = Anum[j] / (Aden[j]  * 1.0);
		BLAS_cp(bm, Y2, i_one, Z2, i_one);
		BLAS_axpy(bm, factor, X2, i_one, Z2, i_one);
		X2 += m;
		Y2 += m;
		Z2 += m;

		
		/* update of r */
		factor = alpha * factor;
		BLAS_cp(bm, Y1, i_one, Z1, i_one);
		BLAS_axpy(bm, factor, X1, i_one, Z1, i_one);
		X1 += m;
		Y1 += m;
		Z1 += m;
	}
#endif

#if 1
	/* logging */
	if ( log && sync->log ) {
		//log_record(sync, dotid, EVENT_ASYNC_FRACTION, pcnt, 0.0);
	}
#endif
}


void __t_scal_cpaxpy_comb_async_release(int dotid, int l, async_t *sync, int p, int bm, int bn, int m, int n, 
	fp_t alpha, fp_t *Anum, fp_t *Aden, fp_t *X1, fp_t *X2, fp_t *Y1, fp_t *Y2, fp_t *Z1, fp_t *Z2) 
{
	int pcnt; 
	pthread_mutex_t *mutex = &sync->mutex;
	int log = 0;

	BLAS_cp(bm, Y2, i_one, Z2, i_one);
	fp_t len = BLAS_dot(bm, X2, I_ONE, Y1, I_ONE);

	pthread_mutex_lock(mutex);
	if ( sync->create != dotid || sync->dot_control != dotid ) {
		++sync->wait;		
		pthread_cond_wait(&sync->cond, mutex);
	} else if ( sync->consume != dotid ) {
		sync->consume = dotid;

		pcnt = sync->pcnt;
		log = 1;
	}
	pthread_mutex_unlock(mutex);


#if 1
	int j;
	for ( j=0; j<bn; ++j) {
		/* update of x */
		fp_t factor = Anum[j] / (Aden[j]  * 1.0);
		BLAS_cp(bm, Y2, i_one, Z2, i_one);
		BLAS_axpy(bm, factor, X2, i_one, Z2, i_one);
		X2 += m;
		Y2 += m;
		Z2 += m;

		
		/* update of r */
		factor = alpha * factor;
		BLAS_cp(bm, Y1, i_one, Z1, i_one);
		BLAS_axpy(bm, factor, X1, i_one, Z1, i_one);
		X1 += m;
		Y1 += m;
		Z1 += m;
	}
#endif

#if 1
	/* logging */
	if ( log && sync->log ) {
		//log_record(sync, dotid, EVENT_ASYNC_FRACTION, pcnt, 0.0);
	}
#endif
}


void __t_cpaxpy_comb4_async(int dotid, async_t *sync, int p, int bm, int bn, int m, int n, fp_t *gamma1, fp_t *gamma2, fp_t *delta, fp_t *sigma2, fp_t *P2, fp_t *V2, fp_t *X2, fp_t *R2, fp_t *S,\
	fp_t *sigma1, fp_t *P1, fp_t *V1, fp_t *X1, fp_t *R1) {
	int pcnt; 
	pthread_mutex_t *mutex = &sync->mutex;
	int log = 0;

	pthread_mutex_lock(mutex);
	pcnt = sync->pcnt;
	if ( sync->dot_control == dotid ) {
		++sync->wait;		
		pthread_cond_wait(&sync->cond, mutex);
	} else if ( sync->consume != dotid ) {
		sync->consume = dotid;
		pcnt = sync->pcnt;
	}
	pthread_mutex_unlock(mutex);

	int j;
	for ( j=0; j<bn; ++j) {
		fp_t beta = gamma1[j] / gamma2[j];
		fp_t lsigma1 = sigma1[j] = delta[j] - beta * beta * sigma2[j];

		/* update of P1 */
		BLAS_cp(bm, R2, i_one, P1, i_one);
		BLAS_axpy(bm, beta, P2, i_one, P1, i_one);

		/* update of V1 */
		BLAS_cp(bm, S, i_one, V1, i_one);
		BLAS_axpy(bm, beta, V2, i_one, V1, i_one);
		
		fp_t alpha = gamma1[j] / lsigma1;
		//printf("cpaxpy_comb4: gamma %.16e %.16e sigma %.16e beta %.16e alpha %.16e\n", gamma1[j], gamma2[j], lsigma1, beta, alpha);
		/* update of X1 */
		BLAS_cp(bm, X2, i_one, X1, i_one);
		BLAS_axpy(bm, alpha, P1, i_one, X1, i_one);

		/* update of R1 */
		alpha = - alpha;
		BLAS_cp(bm, R2, i_one, R1, i_one);
		BLAS_axpy(bm, alpha, V1, i_one, R1, i_one);

		P2 += m;
		V2 += m;
		X2 += m;
		R2 += m;
		S += m;
		P1 += m;
		V1 += m;
		X1 += m;
		R1 += m;
	}

	/* logging */
	if ( log && sync->log ) {
		//log_record(sync, dotid, EVENT_ASYNC_FRACTION, pcnt, 0.0);
	}
}


void __t_scal_cpaxpy_comb_async_concurrent(int dotid, int l, async_t *sync, int p, int bm, int bn, int m, int n, 
	fp_t alpha, fp_t *Anum, fp_t *Aden, fp_t *X1, fp_t *X2, fp_t *Y1, fp_t *Y2, fp_t *Z1, fp_t *Z2) 
{
	int j;
	for ( j=0; j<bn; ++j) {
		/* update of x */
		fp_t factor = Anum[j] / (Aden[j]  * 1.0);
//		printf("alpha1: %e, alpha2: %e, factor: %e\n", Anum[j], Aden[j], factor);
		BLAS_cp(bm, Y2, i_one, Z2, i_one);
		BLAS_axpy(bm, factor, X2, i_one, Z2, i_one);
		X2 += m;
		Y2 += m;
		Z2 += m;

		
		/* update of r */
		factor = alpha * factor;
		BLAS_cp(bm, Y1, i_one, Z1, i_one);
		BLAS_axpy(bm, factor, X1, i_one, Z1, i_one);
		X1 += m;
		Y1 += m;
		Z1 += m;
	}
}


