#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "cg_main.h"
//#include <mkl_blas.h>
//#include <mkl_spblas.h>
#include "mkl.h"
#include "vector.h"
#include "sparsematrix.h"
#include "finout.h"


long cg(SparseMatrix *A, double *rhs, double *x, int s, int b, int maxit) {
	int IONE = 1; 
	double DONE = 1.0, DMONE = -1.0;

	int n = A->dim1; 
	double *mat = A->vval;
	int *vptr = A->vptr;
	int *vpos = A->vpos;

	double umbral = 1.0e-8;

	double *z = (double*) malloc(sizeof(double) * n);
	if ( z == NULL ) {
		fprintf(stderr, "error: could not allocate z\n");
		return 1;
	}
	double *res = (double*) malloc(sizeof(double) * n);
	if ( res == NULL ) {
		fprintf(stderr, "error: could not allocate res\n");
		return 1;
	}
	double *p = (double*) malloc(sizeof(double) * n);
	if ( p == NULL ) {
		fprintf(stderr, "error: could not allocate p\n");
		return 1;
	}

	struct timeval start;
	gettimeofday(&start,NULL);
	
	mkl_dcsrsymv ("U", &n, mat, vptr, vpos, x, z);  // z = A * x
	dcopy (&n, rhs, &IONE, res, &IONE);                            // res = b
	daxpy (&n, &DMONE, z, &IONE, res, &IONE);                      // res -= z
	dcopy (&n, res, &IONE, p, &IONE);                              // p = res
	double gamma = ddot (&n, res, &IONE, res, &IONE);                      // gamma = res' * res
	double tol = sqrt (gamma);                                             // tol = sqrt(gamma) = norm (res)
	int iter = 0;
	while ((iter < maxit) && (tol > umbral)) {
		//print_dvector ("seq.log", "p", n, p);
		mkl_dcsrsymv ("U", &n, mat, vptr, vpos, p, z); 	// z = A * p
		//print_dvector ("seq.log", "z", n, z);
		double alpha = gamma / ddot (&n, p, &IONE, z, &IONE);           // alpha = gamma / (p' * z)
		daxpy (&n, &alpha, p, &IONE, x, &IONE);                     // x += alpha * p;
		double rho = alpha;
		alpha = -alpha;
		daxpy (&n, &alpha, z, &IONE, res, &IONE);                     // res -= alpha * z

		alpha = gamma;                                               // recycle alpha ... alpha = gamma
		gamma = ddot (&n, res, &IONE, res, &IONE);                   // gamma = res' * res
		double beta = gamma / alpha;                                       // alpha = gamma / alpha
		dscal (&n, &beta, p, &IONE);                               // p = beta * p
		daxpy (&n, &DONE, res, &IONE, p, &IONE);                    // p += res
		tol = sqrt (gamma);                                          // tol = sqrt(gamma) = norm (res)
		iter++;

		printf ("iter %i tol %20.16e alpha %20.16e beta %20.16e\n",iter, tol, rho, gamma);
	}

	struct timeval stop;
	gettimeofday(&stop,NULL);

	printf ("it %i res %e\n", iter, tol);
	print_dvector ("x.seq", "x", n, x); 

	unsigned long elapsed = (stop.tv_sec - start.tv_sec) * 1000000;
	elapsed += (stop.tv_usec - start.tv_usec);

	free(z); 
	free(res);
	free(p);

	return elapsed;
} 
