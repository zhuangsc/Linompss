#include "matmul_setup.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "genmat.h"
#include "matfprint.h"
#include "fptype.h"
#include "hb.h"
#include "hbext.h"
#include "iohb.h"
#include "hbconvrt.h"
#include "selfsched.h"


extern selfsched_t *sched;
extern int schedc;

int matmul_setup(const char *fname, int m, int n, int k, int b, int d, int c, void **A, void **B, void **C) {
	hbmat_t *lA = malloc(sizeof(hbmat_t));
	hb_reset(lA);

	printf("warn: A=\"%s\" considered CSR, B and C column major\n", fname);

    readHB_newmat(fname, &(lA->m), &(lA->n), &(lA->elemc), &(lA->vptr), &(lA->vpos), (fp_t **)&(lA->vval));
	hb_flag(lA, FRM_CSR);

	*A = hb2hbh_symcsr(lA, b, NULL, 0);

	/* not set correctly in config */
	m = lA->m;
	k = lA->n ;
	/* not yet supported */
	d = b;
	c = b;

	*B = GENMAT(k, n, k);
	*C = malloc(m * n * sizeof(fp_t));

	int kc = ( k + d - 1 ) / d;
	int lc = schedc = ( n + c - 1 ) / c;
	sched = malloc(lc * sizeof(selfsched_t));
	int w;
	for ( w=0; w<lc; ++w ) {  	
		int nok = sched_ipmalloc(-1, kc, &(sched[w]));
		if ( nok ) {
			fprintf(stderr, "err: could not allocate selfsched\n");
			return 1;
		}
	}

	return 0;
}


int matmul_shutdown(int m, int n, int k, int b, int d, int c, void *A, void *B, void *C) {
	hb_free(A);
	free(B);
	free(C);

	return 0;
}
