#include "matfprint.h"


#include <stdio.h>
#include <string.h>
#include <math.h>

#include "fptype.h"


#ifdef SINGLE_PRECISION

#define fprint_csr2mm			fprint_scsr2mm
#define fprint_csc2mm			fprint_scsc2mm
#define fprint_dense2mm			fprint_sdense2mm
#define print_dense2mm			print_sdense2mm
#define fprint_suff_dense2mm	fprint_suff_sdense2mm
#define print_hb				print_shb

#else

#define fprint_csr2mm			fprint_dcsr2mm
#define fprint_csc2mm			fprint_dcsc2mm
#define fprint_dense2mm			fprint_ddense2mm
#define print_dense2mm			print_ddense2mm
#define fprint_suff_dense2mm	fprint_suff_ddense2mm
#define print_hb				print_dhb

#endif


void fprint_csr2mm(const char *fname, int m, const int *vptr, const int *vpos, const fp_t *vval) {
	FILE *f = fopen(fname, "w");
	int offs = vptr[0] == 0 ? 0 : 1;
	int roffs = 1 - offs;
	
	int i;
	for ( i=0; i<m; ++i ) {
		int rowb = vptr[i] - offs;
		int rowe = vptr[i+1] - offs;
		int j;
		for ( j=rowb; j<rowe; ++j ) {
			if ( vval == NULL ) {
				fprintf(f, "%i\t%i\t1\n", i+roffs, vpos[j] + roffs);
			} else {
				fprintf(f, "%i\t%i\t%.16e\n", i+roffs, vpos[j] + roffs, vval[j]);
			}
		}
	}

	fclose(f);
}


/* can be done easier if matrix is symmetric */
void fprint_csc2mm(const char *fname, int m, const int *vptr, const int *vpos, const fp_t *vval) {
	FILE *f = fopen(fname, "w");
	int offs = vptr[0] == 0 ? 0 : 1;
	int roffs = 1 - offs;

	int *front = malloc(m * sizeof(int));
	memcpy(front, vptr, m * sizeof(int));
	
	printf("offs %i roffs %i\n", offs, roffs);
	int done = 0;
	int i = 0;
	while ( done != m ) {
		int j;
		for ( j=0; j<m; ++j ) {
			if ( front[j] >= 0 ) {
				if ( front[j] == vptr[j+1] ) {
					front[j] = -1;
					++done;
				} else {
					int x = front[j] - offs;
					int r = vpos[x] - offs;
					printf("x %i r %i %.16e\n", x, r, vval[x]);
					if ( r == i ) {
						fprintf(f, "%i\t %i\t %.16e\n", r+1, j+1, vval[x]);
						++front[j];  
					}
				}
			}
		}
		++i;
	}

	free(front);

	fclose(f);
}


void fprint_dense2mm(const char *fname, const char *name, int m, int n, const fp_t *A, int lda) {
	FILE *f = fopen(fname, "w");
	if ( f == NULL ) {
		fprintf(stderr, "err: cannot open %s for writing\n", fname);
	}

	print_dense2mm(f, name, m, n, A, lda);

	fclose(f);
}

// column-major
void print_dense2mm(FILE *f, const char *name, int m, int n, const fp_t *A, int lda) {
	printf("warning: writing obj %s\n", name);

	fprintf(f, "# name: %s\n", name);
	fprintf(f, "# type: matrix\n");
	fprintf(f, "# rows: %i\n", m);
	fprintf(f, "# columns: %i\n", n);

	int i;
	for ( i=0; i<m; ++i ) {
		int j;
		for ( j=0; j<n; ++j ) {
				fprintf(f, "%.16e \n", A[j*lda+i]);
//				fprintf(f, "\n");
		}
//		fprintf(f, "\n");
	}
}


void fprint_suff_dense2mm(const char *fname, int suff, const char *name, int m, int n, const fp_t *A, int lda) {
		char buf1[128];
		char buf2[128];
		sprintf(buf1, "%s_%i.mm", name, suff);
		sprintf(buf2, "%s_%i", name, suff);

		fprint_dense2mm(buf1, buf2, m, n, A, lda);
}


/* Don t use this, use print_csc2mm */
#if 0
void print_csc(FILE *f, hbmat_t *A, int offs) {
	int roffs = 1 - offs;
	int m = A->m;
	int n = A->n;
	int elemc = A->elemc;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	fp_t *vval = A->vval;

	int nanflag = 0;
	int j;
	for ( j=0; j<n; ++j ) {
		int b = vptr[j] - offs; 
		int e = vptr[j+1] - offs;
        int i;
        for ( i=b; i<e; ++i ) {
				fp_t val = vval[i];
				nanflag |= isnan(val);
				
				fprintf(f, "%i %i %.16f\n", vpos[i]+roffs, j+1, vval[i]);
        }
	}

	if ( nanflag ) {
		fprintf(stderr, "warn: found NaN\n");
	}
}


void fprint_csc(const char *fname, hbmat_t *A, int offs) {
	FILE *f = fopen(fname, "w");
	if ( f == NULL ) {
		fprintf(stderr, "error: cannot open %s for writing\n", fname);
	}

	print_csc(f, A, offs);

	fclose(f);
}
#endif


void print_hb(FILE *f, const char *name, hbmat_t *A, int full) {
	fprintf(f, "HB %s @ %p\n", name, A);

	int m = A->m;
	int n = A->n;
	int b = A->b;
	fprintf(f, "m %i n %i elemc %i fill %.4f\n", m, n, A->elemc, (fp_t) 100 * A->elemc / (m * n));
	
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	fp_t *vval = A->vval;
	int csr = IS_CSR(A->type);
	int dim = csr ? m : n;

	if ( vptr != NULL ) {
		int j;
		for ( j=0; j<=dim; ++j ) {
			fprintf(f, "%i ", vptr[j]);
		}
	} else {
		printf("NULL ");
	}
	printf("*\n");
	fflush(0);

	int *vdiag = A->vdiag;
	if ( vdiag != NULL ) {
		int j;
		for ( j=0; j<dim; ++j ) {
			fprintf(f, "%i ", vdiag[j]);
		}
		fprintf(f, "+\n");
	}

	if ( full ) {
		int offs = vptr[0] == 0 ? 0 : 1;
		int j;
		for ( j=0; j<dim; ++j ) {
			int start = vptr[j];
			int elemc = vptr[j+1] - start;
			start -= offs; 

			int i;
			for ( i=0; i< elemc; i++ ) {
				printf("%i ", vpos[start+i] );
			}
			printf(" )\n");
		}

		if ( b == 0 ) {
			int j;
			for ( j=0; j<dim; j++ ) {
				int start = vptr[j];
				int elemc = vptr[j+1] - start;
				start -= offs; 

				int i;
				for ( i=0; i< elemc; i++ ) {
					printf("%.8e ", vval[start+i] );
				}
				printf("\n");
			}
		}
	}
}
