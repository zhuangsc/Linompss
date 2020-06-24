#include <stdio.h>
#include "NoFLA_Compute_td_QR_var31a.h"
#include "NoFLA_QR_td_unb_var1.h"
#include "NoFLA_Apply_td_BlockH_var1.h"

#define min( a, b ) ( (a) < (b) ? (a) : (b) )


/* ========================================================================= */
int NoFLA_Compute_td_QR_var31a( int nb_alg, int m_U, int n_U, int m_D, 
        double * buff_U, int ldim_U,
        double * buff_D, int ldim_D,
        double * buff_t,
        double * buff_S, int ldim_S ) {
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
  double  * buff_U11, * buff_U12, * buff_D1, * buff_D2, * buff_t1, * buff_S1,
          * buff_W12;
  int     k, b, mn_U, n_D2;

  /* Some initializations. */
  mn_U = min( m_U, n_U );

  /* Main loop. */
  for( k = 0; k < mn_U; k += nb_alg ) {
    b    = min( nb_alg, mn_U - k );
    n_D2 = n_U - k - b;

    buff_U11 = &( buff_U[ ldim_U * k + k ] );           // U( k, k )
    buff_U12 = &( buff_U[ ldim_U * ( k + b ) + k ] );   // U( k, k+b )
    buff_D1  = &( buff_D[ ldim_D * k + 0 ] );           // D( 0, k )
    buff_D2  = &( buff_D[ ldim_D * ( k + b ) + 0 ] );   // D( 0, k+b )
    buff_t1  = &( buff_t[ k ] );                        // t( k )
    buff_S1  = &( buff_S[ ldim_S * k + 0 ] );           // S( 0, k )
    buff_W12 = &( buff_S[ ldim_S * ( k + b ) + 0 ] );   // S( 0, k+b )

    /* ( [ U11; D1 ], t1, S1 ) := QR_td( [ U11; D1 ] ). */
    NoFLA_QR_td_unb_var1( b, b, m_D,
        buff_U11, ldim_U,
        buff_D1,  ldim_D,
        buff_t1,
        buff_S1,  ldim_S );

    /* [ U12; D2 ] = Apply_td_BlockH( [ U11; D1 ], S1, W12, [ U12; D2 ] ). */
    NoFLA_Apply_td_BlockH_var1( m_D, n_D2, b,
        buff_D1,  ldim_D,
        buff_S1,  ldim_S,
        buff_W12, ldim_S,
        buff_U12, ldim_U,
        buff_D2,  ldim_D );
  }

  return 0;
}

