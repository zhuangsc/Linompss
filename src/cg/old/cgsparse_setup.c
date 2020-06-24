
double *x = (double* ) malloc(sizeof(double) * n);
	if ( x == NULL ) {
		fprintf(stderr, "err: could not allocate x\n");
		return 1;
	}
	InitDoubles (x, n, 1.0, 0.0);
	double *rhs = (double*) malloc(sizeof(double) * n);;
	if ( rhs == NULL ) {
		fprintf(stderr, "err: could not allocate rhs\n");
		return 1;
	}
	InitDoubles (rhs, n, 0.0, 0.0);

	CopyShiftInts (mat.vpos, mat.vpos, mat.vptr[n], 1);
	CopyShiftInts (mat.vptr, mat.vptr, n+1, 1);
	mkl_dcsrsymv ("U", &n, mat.vval, mat.vptr, mat.vpos, x, rhs);
	printf("x %f %f %f %f\n", x[0], x[1], x[2], x[3]);
	sparse_print_CSC(argv[1], &mat, rhs);

	InitDoubles (x, n, 0.0, 0.0);



