#include <stdio.h>
#include <string.h>
#include "qrca_kernels.h"
#include "qrca_utils.h"

#include "NoFLA_QR_td_unb_var1.h"
#include "NoFLA_Apply_td_BlockH_var1.h"


/* Compute QR factorization of [ U; D ], where U is upper triangular, and 
   returns factors tau into vector t, and triangular factors into matrix S.
   Arguments:
     nb_alg                   Block size to use in the factorization.
     U  in/out  m_U x n_U     Top part of matrix to be factorized.
                              Upper triangular. Defined by ( buff_U, ldim_U ).
     D  in/out  m_D x n_U     Bottom part of matrix to be factorized.
                              Dense. Defined by ( buff_D, ldim_D ).
     t  out     n_U           Vector with tau scalars.
                              Defined by ( buff_t ).
     S  out     nb_alg x n_U  Matrix with triangular factors S.
                              Defined by ( buff_S, ldim_S ).
*/


#define min( a, b ) ( (a) < (b) ? (a) : (b) )


void NoFLA_Compute_td_QR_var31b( int nb_alg, int m_U, int n_U, int skip, 
        double * buff_U, 
        double * buff_D, 
        double * buff_t,
        double * buff_S) {

	double  * buff_U11, * buff_U12, * buff_D1, * buff_D2, * buff_t1, * buff_S1,
  		* buff_W12;
	int     k, b, mn_U, n_D2;

	// For distributed, local memories ...
 	//memset(buff_S, 0, sizeof(double) * nb_alg * m_U);
 	
	NoFLA_Zero_strict_lower_part( m_U, n_U, buff_D, m_U );

	int m=m_U-skip;
	buff_U+=skip;

	mn_U = min( m_U, n_U );


	for( k = 0; k < mn_U; k += nb_alg ) {
		b    = min( nb_alg, mn_U - k );
		n_D2 = n_U - k - b;

		buff_U11 = &( buff_U[ m_U * k + k ] );           // U( k, k )
		buff_U12 = &( buff_U[ m_U * ( k + b ) + k ] );   // U( k, k+b )
		buff_D1  = &( buff_D[ m_U * k + 0 ] );           // D( 0, k )
		buff_D2  = &( buff_D[ m_U * ( k + b ) + 0 ] );   // D( 0, k+b )
		buff_t1  = &( buff_t[ k ] );                        // t( k )
		buff_S1  = &( buff_S[ nb_alg * k + 0 ] );           // S( 0, k )
		buff_W12 = &( buff_S[ nb_alg * ( k + b ) + 0 ] );   // S( 0, k+b )

		// ( [ U11; D1 ], t1, S1 ) := QR_td( [ U11; D1 ] ). 
		NoFLA_QR_td_unb_var1( m, b, m_U,
			buff_U11, m_U,
			buff_D1,  m_U,
			buff_t1,
			buff_S1,  nb_alg );

		// [ U12; D2 ] = Apply_td_BlockH( [ U11; D1 ], S1, W12, [ U12; D2 ] ). 
		NoFLA_Apply_td_BlockH_var1( m_U, n_D2, b,
			buff_D1,  m_U,
			buff_S1,  nb_alg,
			buff_W12, nb_alg,
			buff_U12, m_U,
			buff_D2,  m_U );
  	}
}
