#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>

#include "hb.h"
#include "vector.h"

int* get_sdpos(hbmat_t * A)
{
	int n = A->n;
	int *vptr = A->vptr;
	int *vpos = A->vpos;

	int *sdpos = (int*) calloc(n, sizeof(int));

	int j;
	for ( j=0; j<n; j++ ) {
		int jj = j + 1;

		int start = vptr[j];
		int elemc = vptr[j+1] - start;
		--start;
		
		int pos = 0;
		int i;
		while ( i<elemc && pos<jj ) {
			pos = vpos[start++];
			++i;
		}

		if ( pos > jj ) {
			sdpos[j] = pos;
		}
	}

	return sdpos;
}

void hb_print_CSC(char *fname, hbmat_t *A)
{
        int m = A->m;
        int n = A->n;
        int elemc = A->elemc;
        int *vptr = A->vptr;
        int *vpos = A->vpos;
        fp_t *vval = A->vval;

        FILE *f = fopen(fname, "w");
	if ( f == NULL ) {
		fprintf(stderr, "error: cannot open %s for writing\n", fname);
	}

        
        fprintf(f, "# name: r\n");
        fprintf(f, "# type: matrix\n");
        fprintf(f, "# rows: %i\n", elemc);
        fprintf(f, "# columns: 1\n");
        int j;
        for ( j=0; j<elemc; j++ ) {
                fprintf(f, " %i\n", vpos[j]);
        }

        fprintf(f, "\n# name: c\n");
        fprintf(f, "# type: matrix\n");
        fprintf(f, "# rows: %i\n", elemc);
        fprintf(f, "# columns: 1\n");
        for ( j=0; j<n; j++ ) {
                int vc = vptr[j+1] - vptr[j]; 
                int jj;
                for ( jj=0; jj<vc; jj++) {
                        fprintf(f, " %i\n", j+1);
                }
        }

        fprintf(f, "\n# name: v\n");
        fprintf(f, "# type: matrix\n");
        fprintf(f, "# rows: %i\n", elemc);
        fprintf(f, "# columns: 1\n");
        for ( j=0; j<elemc; j++ ) {
                fprintf(f, " %.16f\n", vval[j]);
        }

	fclose(f);
}

void hb_print_CSC2(char *fname, hbmat_t *A) 
{
	int n = A->n;
	int m = A->m;
	int elemc = A->elemc;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	fp_t *vval = A->vval;

	FILE *f = fopen(fname, "w");
	if ( f == NULL ) {
		fprintf(stderr, "error: cannot open %s for writing\n", fname);
	}

//	fprintf(f, "# name: %s\n", vname);
//	fprintf(f, "# type: sparse matrix\n");
//	fprintf(f, "# nnz: %d\n", elemc);
//	fprintf(f, "# rows: %d\n# columns: %d\n", m, n);
	int j;
	for ( j=0; j<n; ++j ) {
		int b = vptr[j]; 
		int e = vptr[j+1];
        int i;
        for ( i=b; i<e; ++i ) {
                //fprintf(f, "%i\t%i\t%.16f\n", vpos[i], j+1, vval[i]);
				fprintf(f, "%i %i %.16f\n", vpos[i]+1, j+1, vval[i]);
        }
	}

	fclose(f);
}

void hb_print_dense( FILE* str, char * name, hbmat_t *A, int force )
{
	fprintf(str, "%s = [];\n", name);

	int * vptr = A->vptr;
	int * vpos = A->vpos;
	fp_t * vval = A->vval;
	int m = A->m;
	int n = A->n;

	int j;
	for ( j = 0; j < n; j++ ) {
		int start = vptr[j];
		int elemc = vptr[j+1] - start;
		--start;
	
		int c = 0;
		int i;
		for ( i = 0; i < m ; i++ ) {
			fp_t b = 0.0;
			int pos = start + c;
			if ( vpos[pos] == i+1 && c < elemc ) {
				if ( vval == NULL || force) {
					b = 1.0;
				} else {
					b = vval[pos];
					if ( b == 0 ) { 
						fprintf(stderr, "warning: found 0 in non-zero entry!\n");
					}
				}
				++c;
			}
			fprintf(str, "%s(%d,%d) = %.14f;\n", name, i+1, j+1, b);
		}
	}
}

void hb_print_struc(FILE* f, const char *name, hbmat_t *Ahb) 
{
	int *vptr = Ahb->vptr;
	int *vpos = Ahb->vpos;
	fp_t *vval = Ahb->vval;
	int m = Ahb->m;
	int n = Ahb->n;
	int elemc = Ahb->elemc;
	int b = Ahb->b;

	fprintf(f, "dense %s %ix%i b %i\n", name, m, n, b);

	vector_t *cfront = vector_create();
	vector_t *celemc = vector_create();
	int j;
	for ( j = 0; j < n; j++ ) {
		vel_t vel;
		vel.i = vptr[j] - 1;
		vector_insert( cfront, vel );

		vel.i = 0;
		vector_insert( celemc, vel );
	}
		

	int c = 0;
	int i;
	for ( i = 0; i < m; i++ ) {
		int j;
		for ( j = 0; j < n ; j++ ) {
			int currc = vector_get( celemc, j ).i;
			int colc = vptr[j+1] - vptr[j];

			if ( currc < colc ) {
				int start = vector_get( cfront, j ).i;
				int row = vpos[start];

				if ( i == row - 1 ) {
					//fp_t val = vval[start];

					fprintf(f, "X ");

					vel_t vel;
					vel.i = start+1;
					vector_insertat ( cfront, vel, j );
					vel.i = currc+1;
					vector_insertat ( celemc, vel, j );
				} else fprintf(f, "  ");
			} else fprintf(f, "  ");
		}
		fprintf(f, "\n");
	} 

	vector_free(cfront);
	vector_free(celemc);
}

void hbb_print_dense( FILE* str, char * name, hbmat_t *A ) 
{
	fprintf(str, "%s = [];\n", name);

	int * vptr = A->vptr;
	int * vpos = A->vpos;
	fp_t **vval = A->vval;
	int M = A->m;
	int N = A->n;
	int b = A->b;

	//fprintf(stderr, "warning: hbb_print_dense: printing matrix %ix%i b %i\n", M, N, b);

	fp_t *debug = vval[0];

	int J;
	for ( J = 0; J < N; J++ ) {
		int start = vptr[J];
		int blockc = vptr[J+1] - start;
		--start;
	
		int C = 0;
		int I;
		for ( I = 0; I < M ; I++ ) {
			int pos = start + C;
			if ( vpos[pos] == I+1 && C < blockc ) {
				fp_t *B = vval[pos];

				int j;
				for ( j = 0; j < b; j++ )  {
					int i;
					for ( i = 0 ; i < b; i++ ) {
						int r = I * b + i + 1;
						int c = J * b + j + 1;
			
						fprintf(str, "%s(%d,%d) = %.14f;\n", name, r, c, B[j*b+i]);
					}
				}
				
				++C;
			}
		}
	}
}

void hb_print(FILE *f, const char *name, hbmat_t *A, int full) 
{
	fprintf(f, "hbmat %s %p\n", name, A);

	int M = A->m;
	int N = A->n;
	int b = A->b;
	fprintf(f, "M %i N %i elemc %i fill %.4f\n", M, N, A->elemc, (fp_t) 100 * A->elemc / (M * N));
	
	int *vptr = A->vptr;
	int *vpos = A->vpos;

	int j;
	for ( j=0; j<=N; j++ ) {
		fprintf(f, "%i ", vptr[j]);
	}
	printf("*\n");

	int *vdiag = A->vdiag;
	if ( vdiag != NULL ) {
		int j;
		for ( j=0; j<N; j++ ) {
			fprintf(f, "%i ", vdiag[j]);
		}
		fprintf(f, "+\n");
	}
	
	if ( full ) {
		for ( j=0; j<N; j++ ) {
			int start = vptr[j];
			int elemc = vptr[j+1] - start;
			--start; 

			int i;
			for ( i=0; i< elemc; i++ ) {
				printf("%i ", vpos[start+i] );
			}
			printf(" )\n");
		}

		if ( b == 0 ) {
			fp_t *vval = A->vval;
			for ( j=0; j<N; j++ ) {
				int start = vptr[j];
				int elemc = vptr[j+1] - start;
				--start; 

				int i;
				for ( i=0; i< elemc; i++ ) {
					printf("%.4f ", vval[start+i] );
				}
				printf("\n");
			}
			printf("\n");
		}
	}
}

int hb_diff(hbmat_t *A, hbmat_t *B) 
{
	int m = A->m;
	int n = A->n;

	if ( m != B->m ) {
		return 1;
	}

	if ( n != B->n ) {
		return 2;
	}

	int elemc = A->elemc;
	if ( elemc != B->elemc ) {
		return 3;
	}

	int *vptr = A->vptr;
	int *Bvptr = B->vptr;
	int j;
	for ( j=0; j<=n; j++ ) {
		if ( vptr[j] != Bvptr[j] ) {
			return 4;
		}
	}

	int *vval = A->vval;
	int *Bvval = B->vval;
	int c;
	for ( c=0; c<elemc; c++ ) {
		if ( vval[c] != Bvval[c] ) {
			return 5;
		}
	}

	return 0;
}

hbmat_t* hb_cp(hbmat_t *A) 
{
	int m = A->m;
	int n = A->n;
	int elemc = A->elemc;

	hbmat_t *cp = malloc( sizeof(hbmat_t) );
	cp->m = m;
	cp->n = n;
	cp->elemc = elemc;
	cp->b = -1;

	cp->vptr = malloc( sizeof(int) * n );
	memcpy(cp->vptr, A->vptr, n * sizeof(int));
	cp->vval = malloc( sizeof(fp_t) * elemc );
	memcpy(cp->vval, A->vval, elemc * sizeof(fp_t));
	cp->vpos = malloc( sizeof(int) * elemc );
	memcpy(cp->vpos, A->vpos, elemc * sizeof(fp_t));

	cp->vdiag = NULL;

	return cp;
}

void hb_sanity_check(const char *title, hbmat_t *A, int is_hbh)
{
	int M = A->m; 
	int N = A->n; 
	int elemc = A->elemc;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	hbmat_t **vval_hbh = A->vval;
	fp_t *vval_fp = A->vval;

	FILE *info;
	info = fopen("/dev/null", "w");

//#define DEBUG info
#define DEBUG stderr

	fprintf(DEBUG, "-------------------%s----------------------\n",title);
	fprintf(DEBUG, "\t\tHB Sanity Check\n");
	fprintf(DEBUG, "M: %d, N: %d, elemc: %d\n", A->m, A->n, A->elemc);
	int i;
	for ( i = 0; i <= M; ++i )
		fprintf(DEBUG, "%d ", vptr[i]);
	fprintf(DEBUG, "\n--------------------vptr---------------------------\n");
	for ( i = 0; i < elemc; ++i ) 
		fprintf(DEBUG, "%d ", vpos[i]);
	fprintf(DEBUG, "\n");
	if ( is_hbh ) {
		for ( i = 0; i < elemc; ++i ) {
			fprintf(DEBUG, "vval[%d]: %p ", i, vval_hbh[i]);
			fprintf(DEBUG, "\n");
			hb_sanity_check("vval_hbh", vval_hbh[i], 0);
		}
	} else {
		for ( i = 0; i < elemc; ++i )
			fprintf(DEBUG, "vval[%d]: %e ", i, vval_fp[i]);
	}
	fprintf(DEBUG, "\n----------------Check complete----------------------\n");
	fclose(info);
}

/* Copy basic info from B to A */
void hb_init(hbmat_t *A, hbmat_t *B) 
{
	hb_reset(A);
	int M = B->m;
	int elemc = B->elemc;
	A->m = A->n = M;
	A->elemc = elemc;
	A->vptr = malloc( M * sizeof(int));
	A->vpos = malloc( (elemc + 1) * sizeof(int));
	A->vval = malloc( (elemc + 1) * sizeof(fp_t));
}

/* Copy basic info from B to A */
void hbh_init(hbmat_t *A, hbmat_t *B) 
{
	hb_reset(A);
	int M = B->m;
	int elemc = B->elemc;
	A->m = A->n = M;
	A->elemc = elemc;
	A->vptr = malloc( M * sizeof(int));
	A->vpos = malloc( (elemc + 1) * sizeof(int));
	A->vval = malloc( (elemc + 1) * sizeof(hbmat_t*));
}

/* Copy basic info from B to A */
void hb_init_basic(hbmat_t *A, hbmat_t *B)
{
	hb_reset(A);
	int M = B->m;
	int elemc = B->elemc;
	A->m = A->n = M;
	A->elemc = elemc;
}

void hb_reset(hbmat_t *A)
{
	A->m = 0; A->n =0; A->elemc = 0;
	A->vptr = 0; A->vpos = 0; A->vval = 0;
	A->vdiag = NULL; 
	A->b = 0; 
	A->trans = 0; A->orig = 0; A->hyper = 0;
	A->orig_row = 0; A->orig_col = 0; A->e_tree = 0;
	A->type = 0;
	A->FACT = 0;
}

fp_t* hb_energy_row(hbmat_t *A, int bsze)
{
	int m = A->m;
	int elemc = A->elemc;
	int *vptr = A->vptr; int *vpos = A->vpos; fp_t *vval = A->vval;

	int entries = (m + bsze -1)/bsze;
	fp_t *mat_energy = calloc(entries+1, sizeof(fp_t));
	fp_t max = 0.0;
	fp_t tol = 0.0;

	int idx;
	int i;
	for ( i = 0, idx = 0; i < m; i += bsze, idx++) {
		int bptr = vptr[i];
		int er = i+bsze > m ? m : i+bsze;
		int eptr = vptr[er];
		int j;
		for ( j = bptr; j < eptr; j++ ) {
			mat_energy[idx] += fabs(vval[j]);
		}
		tol += mat_energy[idx];
		if (isless(max, mat_energy[idx])) {
			max = mat_energy[idx];
		}
	}
	mat_energy[idx] = max;

	/* Analysis */
	fp_t avg = tol/idx;
	fp_t var, svar;
	for ( int i = 0; i < idx; i++ ) {
		var = FP_POW((mat_energy[idx]-avg), 2.0);
	}
	var /= idx;
	svar = FP_SQRT(var);
	fprintf(stdout, "Matrix energy tol: %e, avg: %e, var: %e, svar: %e\n", tol, avg, var, svar);

	return mat_energy;
}

/* Construct an array of block diagonal submatrices */
void hb_sym_diag_block(hbmat_t *src_mat, int bsze, hbmat_t *diagb)
{
	/* Assuming CSR */
	int m = src_mat->m;
	/* Number of subblocks */
	int bs = (m+bsze-1)/bsze;
	int *svptr = src_mat->vptr; int *svpos = src_mat->vpos;
	fp_t *svval = src_mat->vval;
	int i;
	/* Loop for generating all the diagonal blocks*/
	for ( i = 0; i < bs; i++ ) {
		hbmat_t *d = &diagb[i];
		int elemc = 0;
		int brow = i*bsze; int erow = brow+bsze;
		erow = erow > m ? m : erow;
		int dim = erow - brow;
		d->m = d->n = dim;
		// Allocate individual HB structures
		// Note that vpos and vval size are over-estimated
		int *vptr = malloc((dim+1) * sizeof(int));
		int esze = (svptr[erow] - svptr[brow]);
		int *vpos = malloc(esze * sizeof(int));
		fp_t *vval = malloc(esze * sizeof(fp_t));
		int idx;
		int row;
		/* Traverse through rows */
		for ( row = brow, idx = 0; row < erow; row++ ,idx++) {
			vptr[idx] = elemc;
			int pos = svptr[row]; int epos = svptr[row+1];
			while ( pos < epos ) {
				int col = svpos[pos];
				/* Only take the lower triangular part of the matrix */
//				if ( col >= row && col < erow ) { //Upper
//				if ( col >= brow && col < row ) { //Lower
				if ( col >= brow && col < erow ) { //Complete
					vpos[elemc] = col - brow;
					vval[elemc] = svval[pos];
					elemc++;
				}
				pos++; 
			}
		}
		vptr[idx] = elemc;
		d->elemc = elemc;
		d->vptr = vptr;
		//FIXME using realloc to reduce memory consumption
		d->vpos = vpos;
		d->vval = vval;
		//TODO Remove verifications
//		hb_sanity_check("A_hb", d, 0);
		assert(idx == dim);
		assert(d->vptr != NULL && d->vpos != NULL && d->vval != NULL);
	}
}
