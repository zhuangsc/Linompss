#include "geqrf_main.h"


#include "fptype.h"
#include "task_geqrf.h"
#include "task_larfb.h"
#include "task_tsqrt.h"
#include "task_tsmqr.h"


#define PRIOR_GEQRF		999
#define PRIOR_LARFB		99
#define PRIOR_TSQRT		99
#define PRIOR_TSMQR		9


void geqrf(int diagl, int tr, int tc, int bs, int mt, int nt, fp_t **A, fp_t **tau, fp_t **T) 
{
	int rk = 0;
	int k;
	for( k = 0; k < diagl; k++ ) {
		int row=rk / tr;
		int rskip=rk % tr;

		TASK_GEQRF( bs, tr, tc, rskip, A[ mt*k+row ], tr, tau[ mt*k+row ], T[ mt*k+row ], bs, PRIOR_GEQRF);

		int j;
		for( j = k+1; j < nt; j++ ) {
			TASK_LARFB( bs, tr, tc, tc, rskip, A[ mt*k+row ], tr, T[ mt*k+row ], bs, A[ mt*j+row ], tr, PRIOR_LARFB);
		}

		int i;
		for( i = row+1; i < mt; i++ ) {
			TASK_TSQRT( bs, tr, tc, rskip, A[ mt*k+row ], /*ts,*/A[ mt*k+i ], /*ts,*/tau[ mt*k+i ], T[ mt*k+i ] /*,bs*/, PRIOR_TSQRT);

			int j;
			for( j = k+1; j < nt; j++ ) {
				TASK_TSMQR( bs, tr, tc, rskip, A[ mt*k+i ], T[ mt*k+i ], A[ mt*j+row ], A[ mt*j+i ], PRIOR_TSMQR);
			}
		}
		rk+=tc;
	}
}


#undef PRIOR_GEQRF		
#undef PRIOR_LARFB	
#undef PRIOR_TSQRT
#undef PRIOR_TSMQR	
