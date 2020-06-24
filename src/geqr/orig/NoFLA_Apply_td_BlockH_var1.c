#include <stdio.h>
#include "NoFLA_Apply_td_BlockH_var1.h"
#include "qrca_utils.h"


/* Apply a block Householder defined by ( I; U2 ) and triangular factor S
   to matrix ( C1; C2 ). W is used as workspace.
     U2  in      m x k   Bottom part of Householder vectors (square).
                         Defined by ( buff_U2, ldim_U2 ).
     S   in      k x k   Triangular factor S.
                         Defined by ( buff_S, ldim_S ).
     W   wk      k x n   Workspace.
                         Defined by ( buff_W, ldim_W ).
     C1  in/out  k x n   Top part of matrix to update.
                         Defined by ( buff_C1, ldim_C1 ).
     C2  in/out  m x n   Bottom part of matrix to update.
                         Defined by ( buff_C2, ldim_C2 ).
*/

int NoFLA_Apply_td_BlockH_var1( int m, int n, int k,
        double * buff_U2, int ldim_U2,
        double * buff_S,  int ldim_S,
        double * buff_W,  int ldim_W,
        double * buff_C1, int ldim_C1,
        double * buff_C2, int ldim_C2 ) 
{
  double   d_one = 1.0, d_m_one = -1.0;

  if( ( k == 0 )||( n == 0 ) ) {
    return 0;
  }

  /* W = triu( S )' * ( C1 + U2' * C2 ); */
  NoFLA_Copy( k, n,
              buff_C1, ldim_C1,
              buff_W, ldim_W );

  dgemm_( "Transpose", "No transpose", & k, & n, & m, 
          & d_one, buff_U2, & ldim_U2,
                   buff_C2, & ldim_C2,
          & d_one, buff_W, & ldim_W );

  dtrmm_( "Left", "Upper", "Transpose", "Non-unit", 
          & k, & n, & d_one,
          buff_S, & ldim_S,
          buff_W, & ldim_W );

  /* C2 = C2 - U2 * W; */
  dgemm_( "No transpose", "No transpose", & m, & n, & k,
          & d_m_one, buff_U2, & ldim_U2,
                     buff_W, & ldim_W,
          & d_one, buff_C2, & ldim_C2 );

  /* C1 = C1 - W; */
  NoFLA_Axpy( k, n, d_m_one, buff_W, ldim_W,
                             buff_C1, ldim_C1 );

  return 0;
}
