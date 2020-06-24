#include <mkl.h>
#include <string.h>
#include <math.h>

#include "cg_kernels_step.h"
#include "sparsematrix.h"
#include "finout.h"


void stask_zred(int i, int b, int n, int nb, double *z2, double *z, double *rhs, double *mu) {
	int ione = 1;
	double done = 1.0;
	double dmone = -1.0;

	int ub = n - i;
	ub = (ub < b) ? ub : b;
	z2 += i;

	dscal(&n, &dmone, z2, &done);

	int j;
	for ( j = 0; j < nb; j++ ) {
		daxpy(&ub, &dmone, z2, &ione, z, &ione);  	

		z2 += n;
	}

	daxpy(&ub, &dmone, rhs, &ione, z, &ione);  	
	dscal(&n, &dmone, z, &done);

	mu += ddot(&ub, z, &ione, z, &ione);	
}


void stask_zred2(int i, int b, int n, int nb, double *z2, double *z, double *mu) {
	int ione = 1;
	double done = 1.0;
	double dmone = -1.0;

	int ub = n - i;
	ub = (ub < b) ? ub : b;
	z2 += i;

	int j;
	for ( j = 0; j < nb; j++ ) {
		daxpy(&ub, &done, z2, &ione, z, &ione);  	

		z2 += n;
	}

	mu += ddot(&ub, z, &ione, z, &ione);	
}
