#include "syrk_setup.h"

#include <stdlib.h>
#include <stdio.h>

#include "genmat.h"
#include "matfprint.h"
#include "fptype.h"

#define MAX(a,b) ((a)>(b)?(a):(b))

extern int lda;
extern int ldc;

int syrk_setup(int n, int k, int b, void **A, void **C, void **Cc) 
{
	*A = GENMAT(lda, MAX(n,k), lda);

	fp_t *lC = *C = malloc(ldc * n * sizeof(fp_t));
	fp_t *lCc = *Cc = malloc(ldc * n * sizeof(fp_t));
	GENMAT_SYM_FULL(n, lC);

	int i;
	for( i = 0; i < n*n; ++i ) {
		lCc[i] = lC[i];
	}

	return 0;
}


int syrk_shutdown(int n, int k, int b, void *A, void *C, void *Cc) 
{
	free(A);
	free(C);
	free(Cc);

	return 0;
}
