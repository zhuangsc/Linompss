#include <stdlib.h>
#include <stdio.h>
#include "NoFLA_QR_tS_unb_var1.h"

/* ========================================================================= */
int NoFLA_QR_tS_unb_var1( int m_A, int n_A, double * buff_A, int ldim_A,
                          double * buff_t, double * buff_S, int ldim_S ) {
  double  * w; 
  int     info;

  /* Create vector w for workspace. */
  w = ( double * ) malloc( sizeof( double ) * n_A );

  /* Factorize A. */

  dgeqr2_( & m_A, & n_A, buff_A, & ldim_A, buff_t, w, & info );
  if( info != 0 ) {
    printf( "ERROR in NoFLA_QR_tS_unbvar1. Info of dgeqr2: %d\n", info );
    return info;
  }

  /* Build S. */
  dlarft_( "Forward", "Columnwise", & m_A, & n_A, buff_A, & ldim_A, buff_t,buff_S, & ldim_S ); 

  /* Remove vector w. */
  free( w );

  return 0;
}

