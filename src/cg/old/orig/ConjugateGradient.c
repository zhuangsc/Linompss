#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mkl_blas.h>
#include <mkl_spblas.h>
#include <ScalarVectors.h>
#include <SparseMatrices.h>

int main (int argc, char *argv[]) {
	int i;
	double norm = 0.0;
	double *vec = NULL, *sol = NULL;
	SparseMatrix mat;

	CreateSparseMatrixHB (argv[1], &mat, 1);

	CreateDoubles (&vec, mat.dim1);
	CreateDoubles (&sol, mat.dim1);
	InitDoubles (vec, mat.dim1, 1.0, 0.0);

	InitDoubles (sol, mat.dim1, 0.0, 0.0);

	CopyShiftInts (mat.vpos, mat.vpos, mat.vptr[mat.dim1], 1);
	CopyShiftInts (mat.vptr, mat.vptr, mat.dim1+1, 1);
	mkl_dcsrsymv ("U", &mat.dim1, mat.vval, mat.vptr, mat.vpos, vec, sol);

	// Begin of CG
	int IZERO = 0, IONE = 1; 
	double DONE = 1.0, DMONE = -1.0, DZERO = 0.0;
	int n, iter, maxiter;
	double beta, tol, rho, alpha, umbral;
	double *res = NULL, *z = NULL, *d = NULL;

	n = mat.dim1; maxiter = n; umbral = 1.0e-8;
	CreateDoubles (&res, n); CreateDoubles (&z, n); CreateDoubles (&d, n);
	if ( argc > 2 ) {
		maxiter = atoi( argv[2] );
	}

	InitDoubles (vec, n, DZERO, DZERO);
	
	iter = 0;
	mkl_dcsrsymv ("U", &n, mat.vval, mat.vptr, mat.vpos, vec, z);  // z = A * x
	dcopy (&n, sol, &IONE, res, &IONE);                            // res = b
	daxpy (&n, &DMONE, z, &IONE, res, &IONE);                      // res -= z
	dcopy (&n, res, &IONE, d, &IONE);                              // d = res
	beta = ddot (&n, res, &IONE, res, &IONE);                      // beta = res' * res
	tol = sqrt (beta);                                             // tol = sqrt(beta) = norm (res)
	while ((iter < maxiter) && (tol > umbral)) {
		mkl_dcsrsymv ("U", &n, mat.vval, mat.vptr, mat.vpos, d, z); // z = A * d
		rho = beta / ddot (&n, d, &IONE, z, &IONE);                 // rho = beta / (d' * z)
		daxpy (&n, &rho, d, &IONE, vec, &IONE);                     // x += rho * d;
		rho = -rho;
		daxpy (&n, &rho, z, &IONE, res, &IONE);                     // res -= rho * z
		alpha = beta;                                               // alpha = beta
		beta = ddot (&n, res, &IONE, res, &IONE);                   // beta = res' * res
		alpha = beta / alpha;                                       // alpha = beta / alpha
		dscal (&n, &alpha, d, &IONE);                               // d = alpha * d
		daxpy (&n, &DONE, res, &IONE, d, &IONE);                    // d += res
		tol = sqrt (beta);                                          // tol = sqrt(beta) = norm (res)
		printf ("%d : %20.16e %20.16e %20.16e)\n", iter, tol, -rho, beta);
		iter++;
	}
	printf ("Fin(%d) --> (%d,%20.16e)\n", n, iter, tol);

	RemoveDoubles (&res); RemoveDoubles (&z); RemoveDoubles (&d);

	// End of CG
	
	RemoveDoubles (&sol);
	RemoveDoubles (&vec);
	RemoveSparseMatrix (&mat);
}
