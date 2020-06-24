#include <stdio.h>
#include "NoFLA_Apply_BlockH_var1.h"
#include "qrca_utils.h"

#define min( a, b ) ( (a) < (b) ? (a) : (b) )
#define max( a, b ) ( (a) > (b) ? (a) : (b) )


/* ========================================================================= */
int NoFLA_Apply_BlockH_var1( int m, int n, int k,
        double * buff_U1, int ldim_U1, 
        double * buff_U2, int ldim_U2, 
        double * buff_S,  int ldim_S, 
        double * buff_W,  int ldim_W, 
        double * buff_C1, int ldim_C1,
        double * buff_C2, int ldim_C2 ) {
/* Apply a block Householder defined by ( U1; U2 ) and triangular factor S
   to matrix ( C1; C2 ). W is used as workspace.
     U1  k x k  Top part of Householder vectors (lower triangular).
     U2  m x k  Bottom part of Householder vectors (square).
     S   k x k  triangular factor S.
     W   k x n  Workspace.
     C1  k x n  Top part of matrix to update.
     C2  m x n  Bottom part of matrix to update. */

  double   d_one = 1.0, d_m_one = -1.0;

  //// printf( "  3. Starting NoFLA_Apply_BlockH_var22a. " );
  //// printf( " dims: %d x %d x %d \n", m, n, k );

  /* Quick return. */
  if( ( k == 0 )||( n == 0 ) ) {
    return 0;
  }

  /* W = triu( S )' * ( U1' * C1 + U2' * C2 ); */

  //// FLA_Copy( C1, W );
  NoFLA_Copy( k, n,
              buff_C1, ldim_C1,
              buff_W, ldim_W );

  //// FLA_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, 
  ////           FLA_TRANSPOSE, FLA_UNIT_DIAG,
  ////           FLA_ONE, U1, W );
  dtrmm_( "Left", "Lower", "Transpose", "Unit", 
          & k, & n, & d_one,
          buff_U1, & ldim_U1,
          buff_W, & ldim_W );

  //// FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE,
  ////           FLA_ONE, U2, C2, FLA_ONE, W );
  dgemm_( "Transpose", "No transpose", & k, & n, & m, 
          & d_one, buff_U2, & ldim_U2,
                   buff_C2, & ldim_C2,
          & d_one, buff_W, & ldim_W );

  //// FLA_Trmm( FLA_LEFT, FLA_UPPER_TRIANGULAR, 
  ////           FLA_TRANSPOSE, FLA_NONUNIT_DIAG, 
  ////           FLA_ONE, S, W );
  dtrmm_( "Left", "Upper", "Transpose", "Non-unit", 
          & k, & n, & d_one,
          buff_S, & ldim_S,
          buff_W, & ldim_W );

  /* C2 = C2 - U2 * W; */

  //// FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
  ////           FLA_MINUS_ONE, U2, W, FLA_ONE, C2 );
  dgemm_( "No transpose", "No transpose", & m, & n, & k,
          & d_m_one, buff_U2, & ldim_U2,
                     buff_W, & ldim_W,
          & d_one, buff_C2, & ldim_C2 );

  /* C1 = C1 - U1 * W; */

  //// FLA_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, 
  ////           FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
  ////           FLA_MINUS_ONE, U1, W );
  dtrmm_( "Left", "Lower", "No transpose", "Unit", 
          & k, & n, & d_m_one,
          buff_U1, & ldim_U1, 
          buff_W, & ldim_W );

  //// FLA_Axpy( FLA_ONE, W, C1 );
  NoFLA_Axpy( k, n, d_one, buff_W, ldim_W,
                           buff_C1, ldim_C1 );

  return 0;
}

