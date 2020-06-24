#include "qrca_kernels.h"

#include "NoFLA_QR_td_unb_var1.h"
#include "NoFLA_Apply_td_BlockH_var1.h"


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


#define min( a, b ) ( (a) < (b) ? (a) : (b) )


void NoFLA_Compute_td_QR_var31a_hier( int nb_alg, int m_U, int n_U, int skip, 
        double * buff_U, 
        double * buff_D, 
        double * buff_t,
        double * buff_S) {

	int m=m_U-skip;
	buff_U+=skip;



#if 0
	NoFLA_QR_td_unb_var1( m, n_U, m_U,
		buff_U, m_U,
		buff_D,  m_U,
		buff_t,
		buff_S,  nb_alg );
#else
	// ( [ U11; D1 ], t1, S1 ) := QR_td( [ U11; D1 ] ). 
	NoFLA_QR_td_unb_var1_hier( m, n_U, m_U,
		buff_U, m_U,
		buff_D,  m_U,
		buff_t,
		buff_S,  nb_alg );
#endif
}
