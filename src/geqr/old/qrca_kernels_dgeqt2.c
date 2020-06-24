#include "qrca_kernels.h"

#include <stdio.h>
#include <stdlib.h>


#define min( a, b ) ( (a) < (b) ? (a) : (b) )

/* Compute QR of dense matrix A and returns factors tau into vector t,
   and triangular factors into matrix S.
   Arguments:
     t                   Block size to use in the factorization.
     A  in/out  m_A x n_A     Matrix to be factorized.
                              Defined by ( A, ldim_A ).
     t  out     n_A           Vector with tau scalars.
                              Defined by ( buff_t ).
     S  out     t x n_A  Matrix with triangular factors S.
                              Defined by ( buff_S, ldim_S ).
*/

void dgeqt2_task( int t, int m_A, int n_A, int skip, 
        double * A, 
        double * T,
        double * S) 
{
	A+=skip;
	int m=m_A-skip;
	
	double *LDWORK=(double*) malloc( sizeof(double) * n_A * t);
	int mn_A = min( m_A, n_A );

	int k;
	for( k = 0; k < mn_A; k += t ) {
		int b = min( t, mn_A - k );

		int m_ABR = m - k;
		int n_W12 = n_A - k - b;

		double *buff_ABR = &( A[ m_A * k + k ] );          
		double *buff_A12 = &( A[ m_A * ( k + b ) + k ] );  
		double *buff_t1  = &( T[ k ] );                   
		double *buff_S1  = &( S[ t * k + 0 ] );          

		int info;
		dgeqr2_( &m_ABR, &b, buff_ABR, &m_A, buff_t1, LDWORK, &info );
		if( info != 0 ) {
			printf( "ERROR in NoFLA_QR_tS_unbvar1. Info of dgeqr2: %d\n", info );
		}

		dlarft_( "Forward", "Columnwise", &m_ABR, &b, buff_ABR, &m_A, buff_t1, buff_S1, &t); 

		dlarfb_("Left","Transpose","Forward","Columnwise",\
			&m_ABR, &n_W12, &b,\
			buff_ABR, &m_A,
			buff_S1, &t,\
			buff_A12, &m_A,\
			LDWORK, &n_A);
	}

	free(LDWORK);
}
