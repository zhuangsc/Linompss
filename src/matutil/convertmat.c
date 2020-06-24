#include "convertmat.h"

#include <stdlib.h>


float *convert_d2s(int m, int n, double *D, int ldimD) {
	float *S = calloc(n * ldimD, sizeof(float));
	float *ret = S;

	int j;
	for ( j=0; j<n; ++j ) {
		int i;
		for ( i=0; i<m; ++i) {
			S[i] = (float) D[i];
		}
		S += ldimD;
		D += ldimD;
	}

	return ret;
}


double *convert_s2d(int m, int n, float *S, int ldimS) {
	double *D = calloc(n * ldimS, sizeof(double));
	double *ret  = D;

	int j;
	for ( j=0; j<n; ++j ) {
		int i;
		for ( i=0; i<m; ++i) {
			D[i] = (double) S[i];
		}
		S += ldimS;
		D += ldimS;
	}

	return ret;
}
