#ifndef __SPDSOLV_H__
#define __SPDSOLV_H__


#include "fptype.h"


int spdsolv_setup(int check, int m, int n, int b, fp_t **A, fp_t **B, fp_t **X); 
void spdsolv_shutdown(fp_t *A, fp_t *B, fp_t *X);


#endif // __SPDSOLV_H__
