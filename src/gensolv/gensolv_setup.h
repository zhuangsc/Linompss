#ifndef __GENSOLV_SETUP_H__
#define __GENSOLV_SETUP_H__


#include "fptype.h"


int gensolv_setup(int check, int m, int n, int b, fp_t **A, fp_t **B, fp_t **Aorig, fp_t **Borig); 
void gensolv_shutdown(fp_t *A, fp_t *B, fp_t *Aorig, fp_t *Borig);

#endif // __gensolv_H__
