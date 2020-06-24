#include "qr_check.h"

#include "fptype.h"
#include "densutil.h"
#include "matfprint.h"
#include "fplapack.h"

#include <math.h>
#include <stdlib.h>
#include <sys/time.h>


int qr_check(int check, int m, int n, int br, int bc, void **Rb, void *A)
{
	if( check > 0 ) {
		//fp_t *R = get_canonical_from_hierarchical( "R", m, m, n, n, br, bc, Rb[0] );
		fp_t *R = CMH2CM(m, n, br, bc, Rb[0], m);

		fp_t *T = (fp_t *) malloc( n * bc * sizeof(fp_t) );
		if ( T == NULL ) { 
			perror("check: allocation error"); 
			return 1;
		}

		fp_t *work = (fp_t *) malloc( n * bc * sizeof(fp_t) );
		if ( work == NULL ) { 
			perror("check: allocation error"); 
			return 2;
		}

        struct timeval start;
        gettimeofday( &start, NULL );

		fp_t norm_diff = 0.0;
		int info;
		LAPACK_geqrt(m, n, bc, A, m, T, bc, work, info );


		struct timeval stop;
        gettimeofday(&stop,NULL);
        unsigned long elapsed = stop.tv_usec - start.tv_usec;
        elapsed += (stop.tv_sec - start.tv_sec ) * 1000000;

		if (info < 0 ) {
			perror("check: GEQRT() error\n");
			norm_diff = 1.0;
		}
		else {
			norm_diff = DMAT_RELERR(OMPSSBLAS_UPPERTRIANG, m, n, R, A);
		  	printf("check: diff norm %e (%.2f ms)\n", norm_diff, elapsed / (fp_t) 1e3 );

		}

		free(work);
		free(R);
		free(T);

		int ret = ( norm_diff > 0.0001 ) ? 1 : 0;

		return ret;
	}

	return 0;
}
