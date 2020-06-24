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


int matmul_setup(const char *fname, int m, int n, int k, int b, int d, int c, void **A, void **B, void **C) {
	hbmat_t *lA = malloc(sizeof(hbmat_t));
	hb_reset(lA);

	printf("warn: A=\"%s\" considered CSR, B and C column major\n", fname);

    readHB_newmat(fname, &(lA->m), &(lA->n), &(lA->elemc), &(lA->vptr), &(lA->vpos), (fp_t **)&(lA->vval));
	//one2zero(lA);
	hb_flag(lA, FRM_CSR);

    *A = hb2hbh(lA, b, MAT_CSR);
	//PRINT_HB(stdout, "Ahbh", *A, 1);

	m = lA->m;
	k = lA->n ;
	*B = GENMAT(k, n, k);
	*C = malloc(m * n * sizeof(fp_t));

	return 0;
}


int matmul_shutdown(int m, int n, int k, int b, int d, int c, void *A, void *B, void *C) {
	hb_free(A);
	free(B);
	free(C);

	return 0;
}
