#include "lu_check.h"

#include "fplapack.h"
#include "densutil.h"
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

extern int *IPIV;
fp_t lu_check(int check, int m, int n, fp_t *A, fp_t *Aorig)
{
	if( check ) {
		int mn = min( m,n );
  
		int *ipiv = malloc( n*sizeof(int) );
		if (ipiv == NULL) { 
			perror("Error allocating work"); 
			exit(1); 
		}

		int IINFO;
		LAPACK_getrf(m, n, Aorig, m, ipiv, &IINFO);

#if 0
		if (check) {
			printf( "...done (info %d )\n", IINFO);
		}
#endif

		fp_t norm_diff = MAT_NORMDIFF('i', m, n, A, m, Aorig, m);
		if(check) {
			printf( "norm of difference: %.9e\n", norm_diff);
		}

		free(ipiv);

		return norm_diff;
	  }

	  return 0;
}
