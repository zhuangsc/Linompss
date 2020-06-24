#include "ompss_dense_chol.h"

#include "chol_llmain.h"
#include "chol_rlmain.h"
#include "chol_nllmain.h"
#include "chol_nrlmain.h"

#include <stdio.h>


/* hypermatrix */
int ompss_dchol_hll(int mt, int b, int t, double **Ah) {
	return dchol_hll(mt, b, t, Ah);
}

int ompss_schol_hll(int mt, int b, int t, float **Ah) {
	return schol_hll(mt, b, t, Ah);
}

int ompss_dchol_hrl(int mt, int b, int t, double **Ah) {
	return dchol_hrl(mt, b, t, Ah);
}

int ompss_schol_hrl(int mt, int b, int t, float **Ah) {
	return schol_hrl(mt, b, t, Ah);
}


/* column-major */
int ompss_dchol_ll(int n, int b, int t, double *A, int lda) {
	return dchol_ll(n, b, t, A, lda);
}

int ompss_schol_ll(int n, int b, int t, float *A, int lda) {
	return schol_ll(n, b, t, A, lda);
}

int ompss_dchol_rl(int n, int b, int t, double *A, int lda) {
	return dchol_rl(n, b, t, A, lda);
}

int ompss_schol_rl(int n, int b, int t, float *A, int lda) {
	return schol_rl(n, b, t, A, lda);
}


/* nested */
int ompss_dchol_nhll(int mt, int b, int t, double **Ah) {
	return dchol_nll(mt, b, t, t, Ah);
}

int ompss_schol_nhll(int mt, int b, int t, float **Ah) {
	return schol_nll(mt, b, t, t, Ah);
}

int ompss_dchol_nhrl(int mt, int b, int t, double **Ah) {
	return dchol_nrl(mt, b, t, t, Ah);
}

int ompss_schol_nhrl(int mt, int b, int t, float **Ah) {
	return schol_nrl(mt, b, t, t, Ah);
}



