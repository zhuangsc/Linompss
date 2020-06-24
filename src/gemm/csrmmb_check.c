#include "matmul_check.h"

#include "fptype.h"
#include "fpblas.h"
#include "densutil.h"
#include "hb.h"
#include "hbconvrt.h"
#include "fpsblas.h"
#include "matfprint.h"

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>


int matmul_check(int check, int b, int d, int c, int m, int n, int k, int transa, int transb, void *Ahbh, void *B, int ldb, void *C, int ldc) {
	if ( check ) {
		hbmat_t *Acm = HBH2HB(Ahbh);
		fp_t *Bcm = B;
		fp_t *Ccm = C;

		m = Acm->m; /* have not been set in the main routine, if using a sparse A */
		k = Acm->n;
		ldb = k;
		ldc = m;

		//FPRINT_DENSE2MM("Cc.mm", "C", m, n, Ccm, m);

		fp_t *Cx = malloc( sizeof(fp_t) * ldc * n );

		struct timeval start;
		gettimeofday( &start, NULL );

		SBLAS_csrmm("N", m, n, k, FP_ONE, "GLNF", Acm->vval, Acm->vpos, Acm->vptr, Acm->vptr+1, Bcm, ldb, FP_NOUGHT, Cx, ldc);

		//FPRINT_DENSE2MM("Cx.mm", "C", m, n, Cx, m);

		struct timeval stop;
		gettimeofday(&stop,NULL);
		unsigned long elapsed = stop.tv_usec - start.tv_usec;
		elapsed += (stop.tv_sec - start.tv_sec ) * 1000000;

	 	fp_t norm = DMAT_RELERR(OMPSSBLAS_NTRIANG, m, n, Ccm, Cx);

		printf("check: diff norm %e (%.2f ms)\n", norm, elapsed / (double) 1e3 );

		free( Cx );
	}

	return 0;
}
