#include "task_tsqrt.h"

#include "fptype.h"
#include "blas_ext.h"
#include "qr_assist.h"


#ifdef SINGLE_PRECISION
#define __t_tsqrt 			task_stsqrt
#else
#define __t_tsqrt 			task_dtsqrt
#endif


int NoFLA_QR_td_unb_var1( int m_U, int n_U, int m_D, double * buff_U, int ldim_U, double * buff_D, int ldim_D, double * buff_t, double * buff_S, int ldim_S); 


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
void __t_tsqrt( int nb_alg, int m_U, int n_U, int skip, 
        fp_t * buff_U, 
        fp_t * buff_D, 
        fp_t * buff_t,
        fp_t * buff_S, int p) {

#ifdef DOUBLE_PRECISION
	int m=m_U-skip;
	buff_U+=skip;

	int mn_U = m_U < n_U ? m_U : n_U;

	int k;
	for( k = 0; k < mn_U; k += nb_alg ) {
		int left = mn_U - k;
		int b    = left < nb_alg ? left : nb_alg; 
		int n_D2 = n_U - k - b;

		double *buff_U11 = &( buff_U[ m_U * k + k ] );           // U( k, k )
		double *buff_U12 = &( buff_U[ m_U * ( k + b ) + k ] );   // U( k, k+b )
		double *buff_D1  = &( buff_D[ m_U * k + 0 ] );           // D( 0, k )
		double *buff_D2  = &( buff_D[ m_U * ( k + b ) + 0 ] );   // D( 0, k+b )
		double *buff_t1  = &( buff_t[ k ] );                        // t( k )
		double *buff_S1  = &( buff_S[ nb_alg * k + 0 ] );           // S( 0, k )
		double *buff_W12 = &( buff_S[ nb_alg * ( k + b ) + 0 ] );   // S( 0, k+b )

		// ( [ U11; D1 ], t1, S1 ) := QR_td( [ U11; D1 ] ). 
		NoFLA_QR_td_unb_var1( m, b, m_U,
			buff_U11, m_U,
			buff_D1,  m_U,
			buff_t1,
			buff_S1,  nb_alg );


		// [ U12; D2 ] = Apply_td_BlockH( [ U11; D1 ], S1, W12, [ U12; D2 ] ). 
		dapply_td_BlockH_var1( m_U, n_D2, b,
			buff_D1,  m_U,
			buff_S1,  nb_alg,
			buff_W12, nb_alg,
			buff_U12, m_U,
			buff_D2,  m_U );
	}
#endif
}


/* Compute QR factorization of [ U; D ], where U is upper triangular, and
   returns factors tau into vector t, and triangular factors into matrix S.
   Arguments:
     U  in/out  m_U x n_U     Top part of matrix to be factorized.
                              Upper triangular. Defined by ( buff_U, ldim_U ).
     D  in/out  m_D x n_U     Bottom part of matrix to be factorized.
                              Dense. Defined by ( buff_D, ldim_D ).
     t  out     n_U           Vector with tau scalars.
                              Defined by ( buff_t ).
     S  out     mn_U x mn_U   Matrix with triangular factors S.
                              Defined by ( buff_S, ldim_S ). 
                              mn_U = min( m_U, n_U ).
*/
#ifdef DOUBLE_PRECISION
int NoFLA_QR_td_unb_var1( int m_U, int n_U, int m_D, 
                          double * buff_U, int ldim_U,
                          double * buff_D, int ldim_D,
                          double * buff_t,
                          double * buff_S, int ldim_S ) {
  int mn_U       = m_U < n_U ? m_U : n_U;
  int m_D_plus_1 = m_D + 1;

  double *buff_w = ( double * ) malloc( n_U * sizeof( double ) );

	int j;
  for( j = 0; j < mn_U; j++ ) {
    int n_wt2 = n_U - j - 1;

    dlarfg_( & m_D_plus_1, & buff_U[ j * ldim_U + j ], 
             & buff_D[ j * ldim_D + 0 ], & I_ONE, & buff_t[ j ] );


    /* Update the rest of U and D by applying the Householder transform:
         wt2  := u12t + d1' * D2;
         u12t := u12t - tau1      * wt2;
         D2   := D2   - tau1 * d1 * wt2; */
    dcopy_( & n_wt2, & buff_U[ (j+1) * ldim_U + j ], & ldim_U,
                     & buff_w[ 0 ], & I_ONE );

    dgemv_( "Transpose", & m_D, & n_wt2, 
            & FP_ONE, & buff_D[ (j+1) * ldim_D + 0 ], & ldim_D,  
                     & buff_D[ j * ldim_D + 0 ], & I_ONE,
            & FP_ONE, & buff_w[ 0 ], & I_ONE );

    buff_t[ j ] = - buff_t[ j ];

    daxpy_( & n_wt2, & buff_t[ j ],
            & buff_w[ 0 ], & I_ONE,
            & buff_U[ (j+1) * ldim_U + j ], & ldim_U );

    dger_( & m_D, & n_wt2, & buff_t[ j ],
           & buff_D[ j * ldim_D + 0 ], & I_ONE,
           & buff_w[ 0 ], & I_ONE,
           & buff_D[ (j+1) * ldim_D + 0 ], & ldim_D );
 
    buff_t[ j ] = - buff_t[ j ];


    /* Update current column of S:
         sigma11 := tau1;
         s01     := - tau1 * triu( S00 ) * D0' * d1;  */
    buff_S[ j * ldim_S + j ] = buff_t[ j ];

    dgemv_( "Transpose", & m_D, & j,
            & FP_ONE, & buff_D[ 0 * ldim_D + 0 ], & ldim_D,  
                     & buff_D[ j * ldim_D + 0 ], & I_ONE,
            & FP_NOUGHT, & buff_S[ j * ldim_S + 0 ], & I_ONE );

    dtrmv_( "Upper", "No transpose", "Non-unit", & j,
            &( buff_S[ 0 * ldim_S + 0 ] ), & ldim_S,
            &( buff_S[ j * ldim_S + 0 ] ), & I_ONE );

    dscal_( & j, &( buff_S[ j * ldim_S + j ] ),
                 &( buff_S[ j * ldim_S + 0 ] ), & I_ONE );

    dscal_( & j, &FP_MONE, &( buff_S[ j * ldim_S + 0 ] ), & I_ONE );
  }

  free( buff_w );

  return 0;
}


#endif

