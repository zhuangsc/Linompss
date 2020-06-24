#include "spdsolv_check.h"

#include <math.h>
#include <stdio.h>


int spdsolv_check(int check, int m, int n, fp_t *Xstar, fp_t *X) 
{
	if ( check == 0 ) {
		return 0;
	}

	fp_t norm2 = MAT_NORMDIFF('i', m, n, Xstar, m, X, m);

	printf("check: residual norm = %.16e\n", norm2);

	int ret = norm2 > 1e-6;

	return ret;
}
