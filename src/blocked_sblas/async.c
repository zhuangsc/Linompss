#include "async.h"

#include "task_log.h"
#include "async_struct.h"
#include "fptype.h"
#include "fpblas.h"


#ifdef SINGLE_PRECISION

#define	__async_break	sasync_break
#define __async_conv	sasync_conv

#else

#define __async_break 	dasync_break
#define __async_conv	dasync_conv

#endif

async_stat_t __async_break(int it, fp_t *p, fp_t *c, volatile int *ready, int force)
{
	if ( force || *ready == it ) {
		fp_t pnorm = FP_SQRT(*p);
		fp_t cnorm = FP_SQRT(*c);
		if ( isgreaterequal(cnorm, pnorm) ) {
			return STAT_BROKEN;
		} else {
			return STAT_SYNC;
		}
	} else {
		return STAT_AHEAD;
	}
}

async_stat_t __async_conv(int it, fp_t res, fp_t *c, fp_t div, async_t *sync, int force)
{
	volatile int *ready = &sync->ready;
	if ( force || *ready == it ) {
		fp_t cnorm = FP_SQRT(*c);
		cnorm /= div;
		log_locknrecord(sync, it, EVENT_RESIDUAL, force, cnorm);
		if ( isgreaterequal(res, cnorm) ) {
			return STAT_CONVERGED;
		} else {
			return STAT_SYNC;
		}
	} else {
		return STAT_AHEAD;
	}
}
