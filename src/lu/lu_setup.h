#ifndef __LU_SETUP_H__
#define __LU_SETUP_H__

#include "fptype.h"

void lu_setup(fp_t **A, fp_t **Aorig, int **IPIV, int m, int n, int check);

void lu_shutdown(fp_t *A, fp_t *Aorig, int *IPIV);

#endif // __LU_SETUP_H__
