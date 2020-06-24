#ifndef __TRIANGSOLV_H__
#define __TRIANGSOLV_H__


#include "fptype.h"


int trsm_setup(int check, int m, int n, int b, int lda, int ldb, fp_t **A, fp_t **B, fp_t **X); 
void trsm_shutdown(fp_t *A, fp_t *B, fp_t *X);


#endif // __TRIANGSOLV_H__
