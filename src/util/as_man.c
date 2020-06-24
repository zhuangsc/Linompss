#include "as_man.h"

#include "fptype.h"
#include "task_log.h"
#include "async_struct.h"

#include <stdio.h>
#include <stdlib.h>


#define CGAS_MAN_BREAK 0


/* returns: -1 (divergence), 0 (completed), 1 (not completed) */
async_stat_t ASMAN_BREAK(int it, int est, asman_t *asman, void **norms, async_t *sync, int force, double tol, fp_t *residuals) 
{
	int istar = asman->istar;
	int itop = asman->itop;
	int iprev = asman->iprev;
	int *istack = asman->istack;

	fp_t **lnorms = (fp_t**) norms;
	volatile int *state = &sync->ready;
	//printf("CGAS: %i it %8i state %i @ %p\n", force, it, *state, state);
	if ( it == *state ) {
		int frst = istar == -1;
		fp_t *alphap1 = lnorms[est];	
		fp_t sr2norm = FP_SQRT(alphap1[0]);
		fp_t *alphap2 = lnorms[istar];	
		fp_t prevsr2norm = FP_SQRT(alphap2[0]);

		*residuals = sr2norm;
//		printf("iter: %d, residual: %e\n", it, sr2norm);
//		printf("CGAS: %8i state %i est %i res %.12e best %i prev %i force %i\n", it, *state, est, sr2norm, istar, iprev, force);
//		log_locknrecord(sync, it, EVENT_RESIDUAL, force, sr2norm);
		//printf("%i %i %i %i\n", istack[0], istack[1], istack[2], istack[3]);

		fp_t ordr;
		fp_t prevordr;
		if ( isgreater(sr2norm, prevsr2norm) && !frst ) {
#if CGAS_MAN_BREAK
			ordr = FP_LOG10(sr2norm);
			prevordr = FP_LOG10(prevsr2norm);
			fp_t ordrdif = ordr - prevordr;

			if ( ordrdif > 1.0 ) {
//				printf("ordrdiff\n");
				return STAT_BROKEN;
			} else 
#endif
			{ /* ordr ~ prevordr */
				if ( iprev != istar ) {
					//printf("i!=istar\n");
					istack[--itop] = iprev;
				}
			}
		} else {
			//printf("improved %i %i \n", itop, istar);
			
			if ( !frst ) {
				istack[--itop] = istar;
			}

			if ( iprev != istar ) {
				//printf("i!=istar\n");
				istack[--itop] = iprev;
			}

			istar = est;

			//TODO change the convergence criteria
			if ( islessequal(sr2norm, tol) ) {
//				printf("converge sr2norm: %e\n", sr2norm);
				return STAT_CONVERGED;
			}
		}

		if ( itop < 0 ) { 
			fprintf(stderr, "top %i\n", itop);
		}

		asman->istar = istar;
		asman->itop = itop;
		asman->iprev = est;
		
		return STAT_SYNC;
	} else {
		if ( force ) {
			fprintf(stderr, "break %i: %i %i\n", it, est, *state);
			fprintf(stderr, "err: expected end %i!=%i\n", it, *state);
			return STAT_ERROR;
		}
	}

	return STAT_AHEAD;
}
