#include "task_tsmqr.h"


#include "qr_assist.h"
#include "fptype.h"


#ifdef SINGLE_PRECISION
#define __t_tsmqr 			task_stsmqr
#else
#define __t_tsmqr 			task_dtsmqr
#endif


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
void __t_tsmqr( int nb_alg, int m_G, int n_G, int skip, 
      fp_t * buff_D, /*int ldim_D,*/
      fp_t * buff_S, /*int ldim_S,*/
      fp_t * buff_F, /*int ldim_F,*/
      fp_t * buff_G/*, int ldim_G*/, int p) 
{
	double *buff_W = ( double * ) malloc( nb_alg * n_G * sizeof( double ) );
	int ldim_W = nb_alg;

	buff_F+=skip;
	int m=m_G-skip;

	int k;
	for( k = 0; k < n_G; k += nb_alg ) {
		int left = m - k;
		int b = left < nb_alg ? left : nb_alg; 

		double *buff_D1  = &( buff_D[ m_G * k + 0 ] );   // D( 0, k )
		double *buff_S1  = &( buff_S[ nb_alg * k + 0 ] );   // S( 0, k )
		double *buff_F1  = &( buff_F[ m_G * 0 + k ] );   // F( k, 0 )
  
		dapply_td_BlockH_var1( m_G, n_G, b,
			buff_D1, m_G,
			buff_S1, nb_alg,
			buff_W,  ldim_W,
			buff_F1, m_G,
			buff_G,  m_G );
	}

	free( buff_W );
}
