#include "task_larfb.h"


#include "fptype.h"

#include <stdlib.h>


#ifdef SINGLE_PRECISION

#define __t_larfb 			task_slarfb

#else

#define __t_larfb			task_dlarfb

#endif



/* Apply several block Householder transformations defined by (Ui,Si)
   to matrix C from left side:
     for i = 1, n/t
       C := ( I - Ui * Si * Ui' ) * C.
   Arguments:
     t                   Block size to use in the factorization.
     U  in      m_U x n     Matrix with the Householder vectors.
                              Defined by ( V, ldim_U ).
     S  in      t x n  Matrix with triangular factors S.
                              Defined by ( T, ldim_S ).
     C  in/out  m_U x n_C     Matrix to be applied the reflectors.
                              Defined by ( buff_C, ldim_C ).
*/
void __t_larfb( int t, int m, int n, int k, int skip, fp_t * V, int ldv, fp_t * T, int ldt, fp_t * C, int ldc, int p) 
{
	C += skip;
	V += skip;
	m -= skip;

  	fp_t *W = ( fp_t * ) malloc( n * k * sizeof( fp_t ) );
	int mn = m < n ? m : n; 


  	int l;
	for( l = 0; l < mn; l += t ) {
		int left = mn - l;
		int b = t < left ? t : left;
		int m_U21 = m - l; 

		fp_t *V11 = &( V[ ldv * l + l] );          // U( l, l )
		fp_t *T1  = &( T[ ldt * l + 0 ] );          // S( 0, l )
		fp_t *buff_C1  = &( C[ ldc * 0 + l ] );          // C( l,   0 )

		dlarfb_("Left","Transpose","Forward","Columnwise",\
			&m_U21, &n, &b,\
			V11, &ldv,\
			T1, &t,\
			buff_C1, &ldc,\
			W, &n);
  }

  free( W );
}

