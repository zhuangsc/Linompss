#include "qr_assist.h"

#include "blas_ext.h"

#ifdef SINGLE_PRECISION
#define __apply_td_BlockH_var1				sapply_td_BlockH_var1
#else
#define __apply_td_BlockH_var1				dapply_td_BlockH_var1
#endif


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
int __apply_td_BlockH_var1( int m, int n, int k,
        fp_t * buff_U2, int ldim_U2,
        fp_t * buff_S,  int ldim_S,
        fp_t * buff_W,  int ldim_W,
        fp_t * buff_C1, int ldim_C1,
        fp_t * buff_C2, int ldim_C2 ) 
{

  if( ( k == 0 )||( n == 0 ) ) {
    return 0;
  }

  /* W = triu( S )' * ( C1 + U2' * C2 ); */
  NoFLA_Copy( k, n,
              buff_C1, ldim_C1,
              buff_W, ldim_W );

  dgemm_( "Transpose", "No transpose", & k, & n, & m, 
          & FP_ONE, buff_U2, & ldim_U2,
                   buff_C2, & ldim_C2,
          & FP_ONE, buff_W, & ldim_W );

  dtrmm_( "Left", "Upper", "Transpose", "Non-unit", 
          & k, & n, & FP_ONE,
          buff_S, & ldim_S,
          buff_W, & ldim_W );

  /* C2 = C2 - U2 * W; */
  dgemm_( "No transpose", "No transpose", & m, & n, & k,
          & FP_MONE, buff_U2, & ldim_U2,
                     buff_W, & ldim_W,
          & FP_ONE, buff_C2, & ldim_C2 );

  /* C1 = C1 - W; */
  NoFLA_Axpy( k, n, FP_MONE, buff_W, ldim_W,
                             buff_C1, ldim_C1 );

  return 0;
}
