#include "hbconvrt.h"
#include "hbext.h"
#include "fpsblas.h"
#include "matfprint.h"
#include <assert.h>


#ifdef SINGLE_PRECISION

#define __hbh2hb						shbh2hb
#define __hb2hbh_block					shb2hbh_block
#define __hb_sym_expand					shb_sym_expand

#else

#define __hbh2hb						dhbh2hb
#define __hb2hbh_block					dhb2hbh_block
#define __hb_sym_expand					dhb_sym_expand

#endif


hbmat_t* hb_transpose(hbmat_t *A)
{
	hbmat_t *B = (hbmat_t*) malloc(sizeof(hbmat_t));
	int m = A->m;
	int n = A->n;
	int elemc = A->elemc;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	fp_t *vval = A->vval;
	int acc = 0;

	B->m = n; B->n = m; B->elemc = elemc;
	vector_t *trans_vptr, *trans_vpos, *trans_vval;
	trans_vptr = vector_create_size(n);
	trans_vpos = vector_create_size(elemc);
	trans_vval = vector_create_size(elemc);
	vector_clear(trans_vptr); 
	vector_clear(trans_vpos);
	vector_clear(trans_vval);
	vel_t vptr_vel, vpos_vel, vval_vel;

	int j;
	for( j = 0; j < m; j++){
		vptr_vel.i = trans_vpos->elemc + 1;
		vector_insert(trans_vptr, vptr_vel);
		int i;
		for( i = 0; i < n; i++){
			int k;
			for( k = vptr[i]; k < vptr[i+1]; k++){
				if (j == vpos[k-1] - 1){
					acc++;
				 	vpos_vel.i = i + 1;
					vector_insert(trans_vpos, vpos_vel);
					vval_vel.d = vval[k-1];
					vector_insert(trans_vval, vval_vel);
					break;
				}
			}
		}
	}
	vptr_vel.i = trans_vpos->elemc + 1;
	vector_insert(trans_vptr, vptr_vel);

	B->vptr = vector2int(trans_vptr);
	B->vpos = vector2int(trans_vpos);
	B->vval = vector2double(trans_vval);

	return B;
}


hbmat_t* __hbh2hb(hbmat_t *A) 
{
	hbmat_t *B = malloc(sizeof(hbmat_t));
	hb_reset(B);

	int M = A->m; 
	int N = A->n;
	int elemc = A->elemc;
	int* vptr = A->vptr;
	int* vpos = A->vpos;
	hbmat_t** vval = A->vval;
	int offs = vptr[0]==0 ? 0 : 1;
	
	vector_t *b_vptr = vector_create(); 
	vector_t *b_vpos = vector_create(); 
	vector_t *b_vval = vector_create();
	vel_t vel;
	int bs = vval[0]->m; //Block size can be determined by the rows of the first sub-matrix
	int col_counter = 0;

	int csr = hb_CSR(A);
	int dim = csr ? M : N;

	int K;
	for( K = 0; K < dim; ++K ) {
		hbmat_t *subm = vval[vptr[K] - offs]; //Fetch the first sub-matrix in this column
		int width = csr ? subm->m : subm->n;

		int k;
		for( k = 0; k < width; ++k ) {
			++col_counter;
			vel.i = b_vpos->elemc + offs;
			vector_insert(b_vptr, vel);

			int L;
			for( L = vptr[K]; L < vptr[K+1]; ++L ) {
				int li = L - offs;
				int jump = (vpos[li] - offs) * bs;

				subm = vval[li];
				fp_t *subvval = (fp_t*) subm->vval;
				int  *subvptr = subm->vptr;
				int  *subvpos = subm->vpos;

				int ll;
				for( ll = subvptr[k]; ll < subvptr[k+1]; ++ll ) {
					int lll = ll - offs;
					fp_t el = subvval[lll];
					if( el != 0.0 ) {
						vel_t vel;
						vel.i = subvpos[lll] + jump;
						vector_insert(b_vpos, vel);
						
						vel.d = el;
						vector_insert(b_vval, vel);
					}
				}
			}
		}
	}

	vel.i = b_vpos->elemc + offs;
	vector_insert(b_vptr, vel);
	
	B->m = B->n = col_counter;
	B->elemc = b_vpos->elemc;
	B->vptr = vector2int(b_vptr);
	B->vpos = vector2int(b_vpos);
	B->vval = vector2double(b_vval);
	B->vdiag = NULL;

	return B;
}


hbmat_t *csc2csr(hbmat_t *A)
{
	int m = A->m; int n = A->n; int elemc = A->elemc;
	int* vptr = A->vptr; int* vpos = A->vpos;
	fp_t* vval = A->vval;
	hbmat_t* B = malloc(sizeof( hbmat_t ));
	vector_t *vptr_b, **vpos_b, **vval_b;
	vptr_b = vector_create();
	vpos_b = malloc( m * sizeof(vector_t*));
	vval_b = malloc( m * sizeof(vector_t*));
	vector_clear(vptr_b);
	vel_t b_vel;

	int i;
	for(i = 0; i < m; i++){
		vpos_b[i] = vector_create();
		vector_clear(vpos_b[i]);
		vval_b[i] = vector_create();
		vector_clear(vval_b[i]);
	}

	int j;
	for (j = 0; j < n; j++){
		int i;
		for(i = vptr[j]; i < vptr[j+1]; i++){
			b_vel.i = j;
			vector_insert(vpos_b[vpos[i]], b_vel);
			b_vel.d = vval[i];
			vector_insert(vval_b[vpos[i]], b_vel);
		}
	}

	//Update vptr_b
	b_vel.i = 0;
	vector_insert(vptr_b, b_vel);
	for(i = 0; i < m; i++){
		b_vel = vector_get(vptr_b,i);
		b_vel.i += vpos_b[i]->elemc;
		vector_insert(vptr_b, b_vel);
	}

	for(i = 1; i < m; i++){
		vpos_b[0] = vector_append(vpos_b[0],vpos_b[i]->elem, vpos_b[i]->elemc);
		vector_free(vpos_b[i]);
		vval_b[0] = vector_append(vval_b[0],vval_b[i]->elem, vval_b[i]->elemc);
		vector_free(vval_b[i]);
	}

	B->m = m; B->n = n; B->elemc = elemc;
	B->vptr = vector2int(vptr_b);
	B->vpos = vector2int(vpos_b[0]);
	B->vval = vector2double(vval_b[0]);
	free(vpos_b); free(vval_b);

	return B;
}

void hb_csrcsc(hbmat_t *Acsr, hbmat_t *Acsc, int csr2csc) 
{
	int job0[6] = {0, 0, 0, 0, 0, 1}; //CSR2CSC
	int job1[6] = {1, 0, 0, 0, 0, 1}; //CSC2CSR
	int M = Acsr->m;
	fp_t *acsr = Acsr->vval;
	int *ja = Acsr->vpos;
	int *ia = Acsr->vptr;
	fp_t *acsc = Acsc->vval;
	int *ja1 = Acsc->vpos;
	int *ia1 = Acsc->vptr;
	int info;

	if ( csr2csc )
		SBLAS_CSRCSC( job0, &M, acsr, ja, ia, acsc, ja1, ia1, &info);
	else
		SBLAS_CSRCSC( job1, &M, acsr, ja, ia, acsc, ja1, ia1, &info);
}

void hbh_csrcsc(hbmat_t *Acsr, hbmat_t *Acsc, int csr2csc) 
{

	int job0[6] = {0, 0, 0, 0, 0, 1}; //CSR2CSC
	int job1[6] = {1, 0, 0, 0, 0, 1}; //CSC2CSR
	int M = Acsr->m;
	int N = Acsr->n;
	int elemc = Acsr->elemc;
	Acsc->m = M; Acsc->n = N; Acsc->elemc = elemc;
	fp_t *acsr = Acsr->vval;
	int *ja = Acsr->vpos;
	int *ia = Acsr->vptr;
	fp_t *acsc = Acsc->vval;
	int *ja1 = Acsc->vpos;
	int *ia1 = Acsc->vptr;
	int info;

	if ( csr2csc )
		SBLAS_CSRCSC( job0, &M, acsr, ja, ia, acsc, ja1, ia1, &info);
	else
		SBLAS_CSRCSC( job1, &M, acsr, ja, ia, acsc, ja1, ia1, &info);
}

void* __hb2hbh_block(int I, int J, hbmat_t *A, int b, hbmat_t *Bp) 
{
	int alloc = Bp == NULL;
	
	if ( b < 0 ) {
		fprintf(stderr, "err: b must be positive\n");
		return NULL;
	}

	int m = A->m; int n = A->n;
	int* vptr = A->vptr; 
	int* vpos = A->vpos; 
	fp_t* vval = A->vval;
	int offs = vptr[0] == 0 ? 0 : 1;
	int csr = hb_CSR(A);

	int brow = I*b;  
	int bcol = J*b;
	int rleft = m - brow;
	int cleft = n - bcol;
	int rows = b > rleft ? rleft : b;
	int cols = b > cleft ? cleft : b;
	int erow = brow + rows;
	int ecol = bcol + cols;
	int dimb = csr ? brow : bcol;
	int dime = csr ? erow : ecol;
	int rngb = csr ? bcol : brow;
	int rnge = csr ? ecol : erow;

	vector_t* ab_vptr = vector_create();
	vector_t* ab_vpos = vector_create();
	vector_t* ab_vval = vector_create();
	vel_t vel;

	int L;
	for ( L = dimb; L < dime; ++L ) {
		vel.i = ab_vpos->elemc + offs; 
		vector_insert(ab_vptr, vel);
	
		int k;
		for ( k = vptr[L]; k < vptr[L+1]; ++k ) {
			int lk = k - offs;
			int c = vpos[lk] - offs; 

			if ( c >= rngb && c < rnge ) {
				vel.i = c - rngb + offs;
				vector_insert(ab_vpos, vel);
				vel.d = vval[lk];
				vector_insert(ab_vval, vel);
			}
		}
	}

	vel.i = ab_vpos->elemc + offs;
	vector_insert(ab_vptr, vel);

	if ( alloc ) {
		Bp = malloc(sizeof(hbmat_t));
		hb_reset(Bp);
	}

	if ( ab_vpos->elemc ) {
		Bp->m = rows; 
		Bp->n = cols; 
		Bp->elemc = ab_vpos->elemc;
		Bp->vdiag = NULL;
		Bp->vptr = vector2int(ab_vptr); 
		Bp->vpos = vector2int(ab_vpos);
		Bp->vval = vector2double(ab_vval);
	} else {
		vector_free(ab_vptr);
		vector_free(ab_vpos);
		vector_free(ab_vval);

		return NULL;
	}

	return Bp;
}


hbmat_t* hb2hbh(hbmat_t *A, int b, int is_csr) 
{
	int m = A->m; 
	int n = A->n; 
	int elemc = A->elemc;
	int *vptr = A->vptr; 
	int *vpos = A->vpos; 
	fp_t* vval = A->vval;
	int M = (m+b-1) / b;
	int N = (n+b-1) / b;
	int num = M * N;
	int offs = vptr[0] == 0 ? 0 : 1;

	hbmat_t* hyper = malloc(sizeof(hbmat_t));
	hb_reset(hyper);
	hyper->m = M; hyper->n = N; hyper->vdiag = NULL;
	hyper->orig = A;
	hyper->vval = malloc(num * sizeof(hbmat_t*));
	hbmat_t** hbmat_array = malloc(num * sizeof(hbmat_t*));

	vector_t* ab_vptr = vector_create(); 
	vector_t* ab_vpos = vector_create();
	vel_t pos_val;

	int acc0 = 0;
	int acc = 0;
	int I, J;
	if ( is_csr ) {
		for ( I = 0; I < M; ++I ) {
			pos_val.i = ab_vpos->elemc + offs;
			vector_insert(ab_vptr, pos_val);
			for ( J = 0; J < N; ++J ) {
				hbmat_t *B = __hb2hbh_block(I, J, A, b, NULL);  
				if ( B != NULL ) {
					pos_val.i = J + offs;
					vector_insert(ab_vpos, pos_val);
					((hbmat_t**)hyper->vval)[acc0] = B;
					++acc0;
				}
				++acc;
			}
		}
	} else {
		printf("warn: hb2hbh for csc not yet implemented\n");
	}

	pos_val.i = ab_vpos->elemc + offs;
	vector_insert(ab_vptr, pos_val);
	hyper->elemc = ab_vpos->elemc;
	hyper->vptr = vector2int(ab_vptr);
	hyper->vpos = vector2int(ab_vpos);

	hb_setdiag(hyper);

	return hyper;
}

#if SINGLE_PRECISION
void hyper_sym_csr_task0(int I, int J, hbmat_t *A, int b, int *etree, int *entry)
{
	int m = A->m;
	int n = A->n;
	int* vptr = A->vptr; 
	int* vpos = A->vpos; 
	int brow = I*b; 
	int erow = (I+1)*b;
	erow = erow > m ? m : erow;
	int bcol = J*b; 
	int ecol = (J+1)*b;
	ecol = ecol > n ? n : ecol;
	int offs = vptr[0] == 0 ? 0 : 1;

	*entry = 0;

	if ( etree != NULL ) {
	/* TODO Add path convergence check if necessary */
		int L;
		for ( L = brow; L < erow; ++L ) {
			int k;
			for ( k = vptr[L] - offs; k < vptr[L+1] - offs; ++k ) {
				int ccol = vpos[k] - offs;
				/* 0-offset indexing for cols and row */
				while ( ccol < ecol && ccol <= L ) {
					if ( ccol >= bcol && ccol < ecol ) {
						*entry = 1;
						return;
					}
					ccol = etree[ccol];
				}
			}
		}
	} else {
		int L;
		for ( L = brow; L < erow; ++L ) {
			int k;
			for ( k = vptr[L] - offs; k < vptr[L+1] - offs; ++k ) {
				int ccol = vpos[k] - offs;
				/* 0-offset indexing for cols and row */
				if ( ccol >= bcol && ccol < ecol ) {
					*entry = 1;
					return;
				}
			}
		}

	}
}
#endif


#define FILLINS 1
hbmat_t* hb2hbh_symcsr(hbmat_t *A, int b, int *etree, int alloc) 
{
	int m = A->m; 
	int n = A->n; 
	int elemc = A->elemc;
	int *vptr = A->vptr; 
	int *vpos = A->vpos; 
	fp_t* vval = A->vval;
	int M = ( m + b - 1 ) / b;
	int N = ( n + b - 1 ) / b;
	int num = M * N;

	hbmat_t* H = malloc(sizeof(hbmat_t));
	hb_reset(H);
	H->m = M; 
	H->n = N; 
	H->b = b;
	H->vval = malloc(num * sizeof(hbmat_t*));
	H->e_tree = etree;
	int* hentry = malloc(num * sizeof(int));
	H->orig = A;

	vector_t* ab_vptr = vector_create(); 
	vector_t* ab_vpos = vector_create();

	int acc = 0;
	int I;
	for ( I = 0; I < M; ++I ) {
		int J;
		for ( J = 0; J < N; ++J) {
			hyper_sym_csr_task0(I, J, A, b, etree, &(hentry[acc]));
			++acc;
		}
	}

	acc = 0;
	int acc0 = 0;
	int i;
	for ( i=0, I=0; i < m; i+=b, ++I ) {
		vel_t vel;
		vel.i = ab_vpos->elemc;
		vector_insert(ab_vptr, vel);

		int left = m - i;
		int bi =  b > left ? left : b;

		hbmat_t **hvval = (hbmat_t**) H->vval;
		int j, J;
		for ( j=0, J=0; j < n; j+=b, ++J ) {
			if ( hentry[acc] ) {
				int left = n - j;
				int bj =  b > left ? left : b;

				vel_t vel;
				vel.i = J;
				vector_insert(ab_vpos, vel);

				hbmat_t *Hb = hvval[acc0] = malloc(sizeof(hbmat_t));
				hb_reset(Hb);

				Hb->m = bi;
				Hb->n = bj;
				Hb->hyper = H;
				Hb->orig = A;
				Hb->orig_row = I;
				Hb->orig_col = J;
				Hb->e_tree = etree;
				++acc0;
			}
			++acc;
		}
	}

	int helemc = ab_vpos->elemc;
	if ( alloc ) {
		int vptr_unit = b + 1;
		int vpos_unit = ceil(b * b * FILLINS);
		int vval_unit = vpos_unit;
		int num_vptr = vptr_unit * helemc;
		int num_vpos = vpos_unit * helemc;
		int num_vval = num_vpos;

		pthread_mutex_t *mutexhb = H->mtx = malloc(sizeof(pthread_mutex_t));
		pthread_mutex_init(mutexhb, NULL);
		H->vptr_unit = vptr_unit;
		H->vpos_unit = vpos_unit;
		H->vval_unit = vval_unit;
		H->vptr_pp = 0;
		H->vpos_pp = 0;
		H->vval_pp = 0;
		H->vptr_pool = malloc(num_vptr * sizeof(int));
		H->vpos_pool = malloc(num_vpos * sizeof(int));
		H->vval_pool = malloc(num_vval * sizeof(fp_t));
	}

	vel_t vel;
	vel.i = helemc;
	vector_insert(ab_vptr, vel);
	H->elemc = helemc;
	H->vptr = vector2int(ab_vptr);
	H->vpos = vector2int(ab_vpos);

	free(hentry);

	return H;
}


/*	Expand an symmetric half matrix B to its full counterpart A */
void __hb_sym_expand(hbmat_t *A, hbmat_t *B)
{
	hb_init_basic(A, B);
	int m = A->m;
	A->elemc = B->elemc * 2 - m;
	int nnz = A->elemc;
	A->vptr = malloc((m+1) * sizeof(int));
	A->vpos = malloc(nnz * sizeof(int));
	A->vval = malloc(nnz * sizeof(fp_t));

	int *vptra = A->vptr; int *vposa = A->vpos; fp_t *vvala = A->vval;
	int *vptrb = B->vptr; int *vposb = B->vpos; fp_t *vvalb = B->vval;

	_sn_t *ll_mat = malloc(m * sizeof(_sn_t));

	int i;
	for ( i = 0; i < m; i++ ) {
		ll_mat[i].row = m;
		ll_mat[i].col = -1;
		ll_mat[i].val = 0.0;
		ll_mat[i].next = NULL;
		ll_mat[i].current = &ll_mat[i];
	}

	int vptr_c = 0;
	int elemc_c = 0;

	for ( i = 0; i < m; i++ ) {
		vptra[i] = vptr_c;
		int bptr = vptrb[i]; 
		int eptr = vptrb[i+1];

		/* Fill the lower csr info */
		_sn_t *c = &ll_mat[i];
		while ( c != NULL && c->col != -1 ) {
			vposa[elemc_c] = c->col;
			vvala[elemc_c] = c->val;
			elemc_c++;
			vptr_c++;
			c = c->next;
		}
					
		/* Copy the upper csr info */
		int j;
		for ( j = bptr; j < eptr; j++ ) {
			int col = vposb[j];
			fp_t val = vvalb[j];
			vposa[elemc_c] = col;
			vvala[elemc_c] = val;
			elemc_c++;
			/*-------------------------------------------------*
			 * linked list insert
			 *-------------------------------------------------*/
			_sn_t *head = &ll_mat[col];
			_sn_t *current = ll_mat[col].current;
			if ( current->col != -1 ) {
				current->next = malloc(sizeof(_sn_t));
				current = current->next;
				head->current = current;
			}
			current->row = col;
			current->col = i;
			current->val = val;
			current->next = NULL;
			current->current = current;
			/*-------------------------------------------------*
			 * linked list insert end
			 *-------------------------------------------------*/
		}
		vptr_c += eptr - bptr;
	}

	vptra[m] = vptr_c;

	for ( i = 0; i < m; i++ ) {
		_sn_t *c= ll_mat[i].next;
		while ( c != NULL ) {
			_sn_t *n = c->next;
			free(c);
			c = n;
		}
	}
	free(ll_mat);
}
