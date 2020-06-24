#include <stdlib.h>
#include <stdio.h>
#include "NoFLA_Apply_td_QT_var31a.h"
#include "NoFLA_Apply_td_BlockH_var1.h"
#include "NoFLA_Apply_td_BlockH_var1_hier.h"


#define min( a, b ) ( (a) < (b) ? (a) : (b) )

/* Apply several block Householder transformations defined by (Di,Si) and 
   stored in ( [ I; D ], S ) to matrix [ F; G ] from left side.
     for i = 1, n_D/nb_alg
       [ F; G ] := ( I - [ I; Di ] * Si * [ I; Di ]' ) * [ F; G ].  
   Arguments:
     nb_alg                   Block size to use in the factorization.
     D  in      m_G x n_D     Bottom part of matrix with Householder vectors.
                              Defined by ( buff_D, ldim_D ).
     S  in      nb_alg x n_D  Matrix with triangular factors S.
                              Defined by ( buff_S, ldim_S ).
     F  in/out  n_D x n_G     Top part of matrix to be updated.
                              Defined by ( buff_F, ldim_F ).
     G  in/out  m_G x n_G     Bottom part of matrix to be updated.
                              Defined by ( buff_G, ldim_G ).
*/

void NoFLA_Apply_td_QT_var31a_hier( int nb_alg, int m_G, int n_G, int skip, 
      double * buff_D, /*int ldim_D,*/
      double * buff_S, /*int ldim_S,*/
      double * buff_F, /*int ldim_F,*/
      double * buff_G/*, int ldim_G*/ ) 
{
	double  *buff_D1, * buff_S1, * buff_F1, * buff_W;
	int      ldim_W;

	buff_W = ( double * ) malloc( nb_alg * n_G * sizeof( double ) );
	ldim_W = nb_alg;

	buff_F+=skip;
	int m=m_G-skip;
	

	int k;
	for( k = 0; k < n_G; k += nb_alg ) {
		int b = min( nb_alg, m - k );

		buff_D1  = &( buff_D[ m_G * k + 0 ] );   // D( 0, k )
		buff_S1  = &( buff_S[ nb_alg * k + 0 ] );   // S( 0, k )
		buff_F1  = &( buff_F[ m_G * 0 + k ] );   // F( k, 0 )
  
		NoFLA_Apply_td_BlockH_var1_hier( m_G, n_G, b,
			buff_D1, m_G,
			buff_S1, nb_alg,
			buff_W,  ldim_W,
			buff_F1, m_G,
			buff_G,  m_G );
	}

	free( buff_W );
}
