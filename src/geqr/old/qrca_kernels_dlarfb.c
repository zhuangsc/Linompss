#include "qrca_kernels.h"
#include "qrca_utils.h"

#include <stdlib.h>


#define min( a, b ) ( (a) < (b) ? (a) : (b) )


/* Apply several block Householder transformations defined by (Ui,Si)
   to matrix C from left side:
     for i = 1, n_U/t
       C := ( I - Ui * Si * Ui' ) * C.
   Arguments:
     t                   Block size to use in the factorization.
     U  in      m_U x n_U     Matrix with the Householder vectors.
                              Defined by ( buff_U, ldim_U ).
     S  in      t x n_U  Matrix with triangular factors S.
                              Defined by ( buff_S, ldim_S ).
     C  in/out  m_U x n_C     Matrix to be applied the reflectors.
                              Defined by ( buff_C, ldim_C ).
*/

void dlarfb_task( int t, int m_U, int n_U, int n_C, int skip, 
        double * buff_U, 
        double * buff_S, 
        double * buff_C) 
{
	buff_C+=skip;
	buff_U+=skip;
	int m=m_U-skip;

	//printf( "  3. NoFLA_Apply_dense_QT_var31a\n" );

  	double *buff_W = ( double * ) malloc( n_U * n_U * sizeof( double ) );
	int mn_U = min( m_U, n_U );


  	int k;
	for( k = 0; k < mn_U; k += t ) {
		int b = min( t, mn_U - k );
		int m_U21 = m - k; 

    		double *buff_U11 = &( buff_U[ m_U * k + k] );          // U( k, k )
		double *buff_S1  = &( buff_S[ t * k + 0 ] );          // S( 0, k )
		double *buff_C1  = &( buff_C[ m_U * 0 + k ] );          // C( k,   0 )

	//print_matrix_hierarchical(stdout, "A", m_U, m_U, n_U, n_U, m_U, n_U, buff_C);

	//	printf("dlarfb %i %i %i\n",m_U21,b,n_U);
		dlarfb_("Left","Transpose","Forward","Columnwise",\
			&m_U21, &n_C, &b,\
			buff_U11, &m_U,\
			buff_S1, &t,\
			buff_C1, &m_U,\
			buff_W, &n_U);

#if 0
		dlarfb_("Left","Transpose","Forward","Columnwise",\
			&m_U21, &n_U, &b,\
			buff_U11, &m_U,
			buff_S1, &t,\
			buff_C1, &m_U,\
			buff_W, &n_U);
#endif

#if 0
    NoFLA_Apply_BlockH_var1( m_U21, n_C, b,
        buff_U11, ldim_U,
        buff_U21, ldim_U,
        buff_S1,  ldim_S,
        buff_W,   ldim_W,
        buff_C1,  ldim_C,
        buff_C2,  ldim_C );
#endif
  }

  free( buff_W );
}

