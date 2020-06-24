#include "tasks_potrf.h"


#include "fptype.h"
#include "fplapack.h"
#include "fpmatr.h"


#ifdef SINGLE_PRECISION		

#define __t_potrf		task_spotrf

#else

#define __t_potrf		task_dpotrf

#endif


void __t_potrf(ompsslapack_t uplo, int n, fp_t *A, int lda, int p) 
{
	int info;
	LAPACK_potrf(uplo, n, A, lda, info);
}
