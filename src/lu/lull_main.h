#ifndef __LULL_MAIN_H__
#define __LULL_MAIN_H__

#include "fptype.h"

#ifdef SINGLE_PRECISION
#define LULL_MAIN		slull_main
#endif

#ifdef DOUBLE_PRECISION
#define LULL_MAIN		dlull_main
#endif

void slull_main(int m, int n, int b, float *A, int lda, int *IPIV);
void dlull_main(int m, int n, int b, double *A, int lda, int *IPIV);


#endif // __LULL_MAIN_H__
