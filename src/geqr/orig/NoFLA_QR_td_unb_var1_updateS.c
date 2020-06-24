#include "qrca_kernels.h"


void NoFLA_QR_td_unb_var1_updateS( int mn_U, int m_D, 
                          double * buff_D, int ldim_D,
                          double * buff_t,
                          double * buff_S, int ldim_S ) {

  	double d_one = 1.0;
  	double d_zero = 0.0; 
	double d_minus_one = -1.0;
  	int i_one = 1;

	/* Update current column of S:
         sigma11 := tau1;
         s01     := - tau1 * triu( S00 ) * D0' * d1;  */
	int j;
	for( j = 0; j < mn_U; j++ ) {
		buff_S[ j * ldim_S + j ] = buff_t[ j ];

		dgemv_( "Transpose", & m_D, & j,
            		& d_one, & buff_D[ 0 * ldim_D + 0 ], & ldim_D,  
                     	& buff_D[ j * ldim_D + 0 ], & i_one,
            		& d_zero, & buff_S[ j * ldim_S + 0 ], & i_one );

    		dtrmv_( "Upper", "No transpose", "Non-unit", & j,
            		&( buff_S[ 0 * ldim_S + 0 ] ), & ldim_S,
            		&( buff_S[ j * ldim_S + 0 ] ), & i_one );

    		dscal_( & j, &( buff_S[ j * ldim_S + j ] ),
	 		&( buff_S[ j * ldim_S + 0 ] ), & i_one );

		dscal_( & j, & d_minus_one, &( buff_S[ j * ldim_S + 0 ] ), & i_one );
  	}
}
