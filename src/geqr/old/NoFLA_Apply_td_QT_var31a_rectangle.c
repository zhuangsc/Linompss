#include <stdlib.h>
#include <stdio.h>
#include "NoFLA_Apply_td_QT_var31a.h"
#include "NoFLA_Apply_td_BlockH_var1.h"

#define min( a, b ) ( (a) < (b) ? (a) : (b) )


void NoFLA_Apply_td_QT_var31a( int nb_alg, int m_G, int n_G, int n_D,
      double * buff_D, int ldim_D,
      double * buff_S, int ldim_S,
      double * buff_F, int ldim_F,
      double * buff_G, int ldim_G ) {
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
  double   * buff_D1, * buff_S1, * buff_F1, * buff_W;
  int      k, b, ldim_W;

  /* Create nb_alg x n_G temporal matrix W. */
  buff_W = ( double * ) malloc( nb_alg * n_G * sizeof( double ) );
  ldim_W = nb_alg;

  for( k = 0; k < n_D; k += nb_alg ) {
    b = min( nb_alg, n_D - k );

    buff_D1  = &( buff_D[ ldim_D * k + 0 ] );   // D( 0, k )
    buff_S1  = &( buff_S[ ldim_S * k + 0 ] );   // S( 0, k )
    buff_F1  = &( buff_F[ ldim_F * 0 + k ] );   // F( k, 0 )
  
    NoFLA_Apply_td_BlockH_var1( m_G, n_G, b,
        buff_D1, ldim_D,
        buff_S1, ldim_S,
        buff_W,  ldim_W,
        buff_F1, ldim_F,
        buff_G,  ldim_G );
  }

  free( buff_W );
}

