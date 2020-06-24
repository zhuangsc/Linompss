#ifndef __CG_FUNCTIONS_STEP_H__
#define __CG_FUNCTIONS_STEP_H__


#include "cg_kernels.h"
#include "cg_kernels_step.h"


static inline __attribute__((always_inline))  c_R_mu(int s, int n, int b, double *mu, double *A, double *R, double *tmP, double *rhs) {
	double *r = &R[0];
	double *p = &R[n];

	// compute r and (r,r)
	int i;
	for ( i = 0; i < n; i+=b ) {
		stask_mv(i, b, n, nb, A, r, &p[i], tmP);
	}
	for ( i = 0; i < n; i+=b ) {
		stask_zred(i, b, n, nb, tmP, &p[i], rhs, mu);
	}

	if ( s > 1 ) {
		r = p;
		p += n;
		// compute A * r and (A * r, A * r)
		int i;
		for ( i = 0; i < n; i+=b ) {
			task_mv(i, b, n, nb, A, r, &p[i], tmP);
		}
		for ( i = 0; i < n; i+=b ) {
			task_zred2(i, b, n, nb, tmP, &p[i], mu);
		}

		// compute A^2 * r etc
		int k;
		for ( k = 2; k <= s; k++ ) {
			double *p = r + n;
			int i;
			for ( i = 0; i < n; i+=b ) {
				task_mv(i, b, n, nb, A, r, &p[i], tmP);
			}
			for ( i = 0; i < n; i+=b ) {
				stask_zred(i, b, n, nb, tmP, &p[i]);
			}
			r += n;
		}
	}
}


#endif // __CG_FUNCTIONS_STEP_H__
