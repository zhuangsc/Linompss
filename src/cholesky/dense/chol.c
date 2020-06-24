#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>


#include "fptype.h"
#include "chol_config.h"
#include "chol_check.h"
#include "chol_llmain.h"
#include "chol_rlmain.h"
#include "chol_nllmain.h"
#include "chol_nrlmain.h"
#include "chol_prlmain.h"
#include "chol_setup.h"

#include "matfprint.h"
#include "file_log.h"


int m; /* order of input matrix */
int ts; /* block size */
int bs; /* sub-block size */
int reps; /* number of repetitions */
int check; /* check result */
int mt; /* number of blocks */
int mr; /* rounded m */
int mtleft;

/* dense */
fp_t *A;
fp_t **Ah;
fp_t *Aorig;	
log_f logf_chol;

int main(int argc, char* argv[]) 
{
	if ( chol_config(argc, argv) ) {
		return 1;
	}
	
	if ( chol_setup(check, m, mr, ts, bs, mt, mtleft) ) {
		fprintf(stderr, "err: allocating matrix\n");
		return 2;
	}

	logf_init(&logf_chol, reps, "chol");
	int r;
	for (r = 0; r < reps; r++) {
		logf_record_stime(&logf_chol);

#if USE_HLL
		CHOL_HLL(mt, ts, bs, Ah);
#elif USE_HRL
		CHOL_HRL(mt, ts, bs, Ah);
#elif USE_LL
		CHOL_LL(m, ts, bs, A, m);
#elif USE_RL
		CHOL_RL(m, ts, bs, A, m);
#elif USE_PRL
		//elapsed += chol_prl(mr, mt, ts, bs, Ah);
#elif USE_NLL
		CHOL_NLL(mt, ts, bs, bs, Ah);
#elif USE_NRL
		CHOL_NRL(mt, ts, bs, bs, Ah);
#endif

		#pragma omp taskwait

		logf_record_etime(&logf_chol);
  	}

 
  	int ret = chol_check(check, m, mr, ts, Ah, A, Aorig);

	logf_dump(&logf_chol);
	chol_shutdown();

	return ret;
}
