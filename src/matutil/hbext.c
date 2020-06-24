#include "hbext.h"

//#include "vector.h"


#if 0
void hb_sync(hbmat_t *A, hbmat_t *B) {
	if ( B->vval == NULL ) {
		fprintf( stderr, "warning: sync allocating vval for B\n");
		size_t s = B->b==0 ? sizeof(fp_t) : sizeof(fp_t*);
		B->vval = (void*) calloc( B->elemc, s );
	}
	fp_t *avval = A->vval;
	fp_t *bvval = B->vval;
	int *avpos = A->vpos;
	int *bvpos = B->vpos;
	int *avptr = A->vptr;
	int *bvptr = B->vptr;


	int m = A->m;
	int n = A->n;

	if ( m!=B->m || n!=B->n ) {
		fprintf(stderr, "hbsync: A and B have different dimensions\n");
		return; 
	}	

	int j;
	for ( j=0; j<n; j++ ) {
		int bstart = bvptr[j] - 1;
		int bend = bvptr[j+1] - 1;
		int belemc = bend - bstart;

		int astart = avptr[j] - 1;
		int aend = avptr[j+1] - 1;
		int aelemc = aend - astart;

		int bcur = 0;
		int i;
		for ( i=0; i<aelemc; i++ ) {
			int arow = avpos[ astart + i ];
			int brow = bvpos[ bstart + bcur];
			while ( brow < arow ) {
				++bcur;
				brow = bvpos[ bstart + bcur];
			}

			bvval[ bstart + bcur ] = avval[ astart + i];
			++bcur;
		}
	}
}
#endif



void hb_markudiag(hbmat_t *A) {
	int elemc = A->elemc;
	int* vptr = A->vptr;
	int *vpos = A->vpos;
	int n = A->n;

	int *udiagc = A->udiagc = (int*) malloc( sizeof(int) * n );

	int j;
	for ( j=0; j<n; j++ ) {
		int start = vptr[j];
		int rowc = vptr[j+1] - start;
		--start;

		int jj = j + 1;
		int i = 0;
		while ( vpos[start+i] < jj ) {
			++i;
		}
		udiagc[j] = i;
	}
}


#if 0
void hbm_setdiag(hbmatm_t *A) {
	int* vptr = A->vptr;
	//dll_t *vpos = A->vpos;
	dll_t *vval = A->vval;
	int n = A->n;

	int *vdiag = A->vdiag = (int*) malloc( sizeof(int) * n );
	dll_t **vdiagl = A->vdiag = (dll_t**) malloc( sizeof(dll_t*) * n );

	int j;
	for ( j=0; j<n; j++ ) {
		int start = vptr[j];  
		int rowc = vptr[j+1] - start;
	
		int c = 0;
		while ( c < rowc ) {
			//vpos = vpos->next;
			vval = vval->next;
			if ( vval->e.i == j+1 ) {
				vdiag[j] = start + c;
				vdiagl[j] = vval;
				++c;
				//printf("diag of col %i located at %i\n", j+1, start+i+1);
				break;
			} //else printf("%i has %i\n", i+1, vpos[start+i]);
			++c;
		}
		while ( c < rowc ) {
			//vpos = vpos->next;
			vval = vval->next;
			++c;
		}
		
	}
}


void csrm_setdiag(hbmatm_t *A) {
#if 0
	int elemc = A->elemc;
	int* vptr = A->vptr;
	int *vpos = A->vpos;
	int m = A->m;
	int n = A->n;

	int *vdiag = A->vdiag = (int*) malloc( sizeof(int) * elemc );

	int i;
	for ( i=0; i<m; i++ ) {
		int start = vptr[i];
		int colc = vptr[i+1] - start;
		--start;
	
		int j;
		for ( j=0; j<colc; j++ ) {
			if ( vpos[start+j] == i+1 ) {
				vdiag[i] = start + j + 1;
	//			printf("diag of col %i located at %i\n", j+1, start+i+1);
				break;
			} //else printf("%i has %i\n", i+1, vpos[start+i]);
		}
	}
#endif
}
#endif

void one2zero(hbmat_t* in_matrix)
{
	int m = in_matrix->m;
	int elemc = in_matrix->elemc;
	int *vptr = in_matrix->vptr;
	int *vpos = in_matrix->vpos;

	int i;
	for ( i = 0; i <= m; i++ ) {
		vptr[i]--;
	}
	for ( i=0; i<elemc; i++ ) {
		vpos[i]--;
	}
}


int* etree(hbmat_t* A)
{
	int n = A->n;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	int *parent = calloc(n, sizeof(int));
	int *ancestor = calloc(n, sizeof(int)); 

	if ( !parent || !ancestor ) {
		return NULL;
	}

	int column;
	for ( column = 0; column < n; column++ ) {
		parent[column] = -1;
		ancestor[column] = -1;

		int p;
		for ( p = vptr[column]; p < vptr[column+1]; p++ ) {
			int inext;
			int i;
			for ( i = vpos[p]; i != -1 && i < column; i = inext ) {
				inext = ancestor[i];
				ancestor[i] = column;
				if (inext == -1)
					parent[i] = column;
			}
		}
	}

	free(ancestor);

	return parent;
}


void hb_setdiag(hbmat_t *A) 
{
	int* vptr = A->vptr;
	int *vpos = A->vpos;
	int m = A->m;
	int n = A->n;
	int csc = IS_CSC(A->type);
	int offs = vptr[0] == 0 ? 0 : 1;

	int *vdiag = A->vdiag = (int*) malloc( sizeof(int) * n );

	int dim1 = csc ? m : n ;
	int dim2 = csc ? n : m ;

	int k;
	for ( k=0; k<dim1; ++k ) {
		int start = vptr[k];
		int c = vptr[k+1] - start;
		start -= offs;;

		vdiag[k] = -1;
	
		int l;
		for ( l=0; l<c; ++l ) {
			if ( vpos[start+l] == k+offs ) {
				vdiag[k] = start + l + offs;
				break;
			}
		}

		if ( vdiag[k] == -1 ) {
			fprintf(stderr, "warn: %p does not have an entry on diag %i\n", A, k+offs);
		}
	}
}
