#include "genmat.h"

#include "fptype.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


// generates an SPD matrix in blocked layout in place
void GENMAT_HYPER_SPD(int n, int tn, int nleft, int b, fp_t *A) {
	fp_t *Ain = &A[0];

	int j;
	for (j = 0; j < tn; j++ ) {
		int jend=j==tn-1;
		int i;
		for( i = 0; i < tn; i++ ) {
			int iend=i==tn-1;
			int jj; 
			for( jj=0 ; jj < b; jj++ ) {
				int ii;
				for ( ii=0 ; ii < b; ii++ ) {
  					int skipright = iend?ii>=b-nleft:0;
  					int skipbotn = jend?jj>=b-nleft:0;
	  
  					if ( !skipright && !skipbotn ) {
    						fp_t dran = drand48();
    						if ( dran == 0.0 )  printf("generated 0\n");
							if( i==j && ii ==jj ) {
								dran += n;
							}
    						*Ain = dran;

    						//printf("A(%i,%i)=%.4f\n",i*b+ii+1,j*b+jj+1,dran);
  					} else {
						if ( i==j & ii==jj ) {
							*Ain=1.00;
						}
					}
					Ain++;
				}
			}
		}
  	}
}

// generates a diagonally dominant matrix in column-major layout
fp_t * GENMAT_SPD(int n, int ldim) 
{
	fp_t *A = malloc(sizeof(fp_t) * ldim * n);
	if ( A == NULL ) { 
		return NULL;
	}

	srand48(time(NULL));

	int j;
	for (j = 0; j < n; ++j ) {
		int i;
		for( i = j; i < n; ++i ) {
			fp_t dran = drand48();
			if ( i==j ) 
				dran += n;
			A[j*ldim+i] = A[i*ldim+j] = dran;
		}
  	}

	return A;
}


// generates a symmetrix matrix in column-major layout 
void GENMAT_SYM_FULL(int m, fp_t *A) 
{
	srand48(time(NULL));

	int j;
	for (j = 0; j < m; ++j ) {
		int i;
		for( i = j; i < m; ++i ) {
			fp_t dran = drand48();
			A[j*m+i] = A[i*m+j] = dran;
		}
  	}
}

// generates a symmetrix matrix in column-major layout 
void GENMAT_SYM_LTRIANG(int m, fp_t *A) 
{
	srand48(time(NULL));
	int j;
	for (j = 0; j < m; ++j ) {
		fp_t dran = drand48();
		A[j*m+j] = dran;
		int i;
		for( i = j+1; i < m; ++i ) {
			fp_t dran = drand48();
			A[j*m+i] = dran;
			A[i*m+j] = 0.0;
		}
  	}
}


fp_t* GENMAT(int m, int n, int ldim) 
{
	if ( ldim < m ) {
		return NULL;
	}

	srand48(time(NULL));
	fp_t *A = malloc(sizeof(fp_t) * ldim * n);
	if ( A == NULL ) { 
		return NULL;
	}

	int j;
	for (j = 0; j < n; ++j ) {
		int i;
		for( i = 0; i < m; ++i ) {
			fp_t dran = (fp_t) drand48();
			A[j*ldim+i] = dran;
		}
  	}

	return A;
}

//In-place random matrix generator column-major layout
void GENMAT_IP(fp_t *A, int m, int n, int scale) 
{
	srand48(time(NULL));

	int j;
	for (j = 0; j < n; ++j ) {
		int i;
		for( i = 0; i < m; ++i ) {
			A[j*m+i] = (fp_t) drand48();
			if ( j == i )
				A[j*m+i] += scale;
		}
  	}
}


void* GENMAT_EYE(int m, int ldim) 
{
	int m2 = m * ldim;
	fp_t *A = malloc(m2 * sizeof(fp_t));

	int j;
	for (j = 0; j < m; ++j ) {
		A[j*ldim+j] = FP_ONE;
  	}

	return A;
}


void* GENMAT_ZERO(int m, int ldim) 
{
	int m2 = m * ldim;
	fp_t *A = calloc(m2, sizeof(fp_t));

	return A;
}



void GENMAT_LTRIANG_ONES(int m, fp_t *A) {
	int m2 = m * m;
	memset(A, 0, m2 * sizeof(fp_t));

	int j;
	for (j = 0; j < m; ++j ) {
		int i;
		for ( i=j; i<m; ++i ) {
			A[j*m+i] = FP_ONE;
		}
  	}
}


// generates a random vector
#if 0
A vector is a matrix with a single column. This is redundant. Use GENMAT instead. 
void GENVEC_RANDOM(fp_t *vector, int dim, fp_t scale) {
	srand48(time(NULL));

	int i;
	for(i = 0; i < dim; ++i){
		vector[i] = drand48() * scale;
	}
}
#endif


fp_t* GENMAT_COSTAS(size_t m, size_t n) 
{
	srand48( ( long int ) time(NULL) );
	//srand48( ( long int ) 666 );
	
	fp_t *A = malloc(m * n * sizeof(fp_t));
	if ( A == NULL ) {
		return NULL;
	}

	fp_t *pA = A;
	size_t s = m * n;
	size_t c;
	for ( c=0; c<s ; ++c ) {
		fp_t a = drand48();
		a = a > 0.5 ? 1.0 : -1.0;
		*pA = a;
		++pA;
	}

	return A;
}


fp_t *GENMAT_SYMDD(size_t m, size_t n) 
{
	srand48( ( long int ) time(NULL) );
	//srand48( ( long int ) 666 );
	
	fp_t *A = malloc(m * n * sizeof(fp_t));
	if ( A == NULL ) {
		return NULL;
	}

	size_t j;
	for (j=0; j<n ; ++j) {
		fp_t *ap = &A[j*m+0];

		size_t i;
		for (i=0; i<j ; ++i) {
			*ap++ = A[i*m+j];
		}    

		fp_t f = drand48() + n;
		*ap++ = f;
		if ( f == 0.0 ) {
			printf("generated 0\n");
		}
		

		for (i=j+1; i<m ; ++i) {
			fp_t f = drand48();
			*ap++ = f;
			if ( f == 0.0 ) {
				printf("generated 0\n");
			}

		}    
	}

	return A;
}


fp_t *GENMAT_COVARIANT(size_t n, fp_t p, fp_t k) 
{
	if ( p <= 0.0 ) {
		fprintf(stderr, "err: a covariant matrix requires p > 0\n");
		return NULL;
	}
	
	fp_t mp = -p;
	fp_t kinv = 1 / k;

	size_t s = n * n;
	fp_t *A = malloc(s * sizeof(fp_t));
	if ( A == NULL ) {
		return NULL;
	}

	size_t j;
	for (j=0; j<n ; ++j) {
		fp_t *ap = &A[j*n+0];

		size_t i;
		for (i=0; i<j ; ++i) {
			fp_t v = FP_POW(j-i, mp);
	//		printf("generated (%i,%i) %e\n", i, j, v);
			*ap++ = v;
		}    

		fp_t v = 1 + FP_POW(j+1, kinv);
	//	printf("generated (%i,%i) %e\n", i, j, v);
		*ap++ = v;

		for (i=j+1; i<n ; ++i) {
			fp_t v = FP_POW(i-j, mp);
	//		printf("generated (%i,%i) %e\n", i, j, v);
			*ap++ = v;
		}    
	}

	return A;
}
