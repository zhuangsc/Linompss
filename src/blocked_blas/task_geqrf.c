#include "task_geqrf.h"


#include "fptype.h"
#include "fplapack.h"

#include <stdio.h>
#include <stdlib.h>


#ifdef SINGLE_PRECISION

#define __t_geqrf 			task_sgeqrf

#else

#define __t_geqrf			task_dgeqrf

#endif


/* Compute QR of dense matrix A and returns factors tau into vector t,
   and triangular factors into matrix S.
   Arguments:
     t                   Block size to use in the factorization.
     A  in/out  m_A x n     Matrix to be factorized.
                              Defined by ( A, ldim_A ).
     tau  out     n           Vector with tau scalars.
                              Defined by ( buff_t ).
     T  out     t x n  Matrix with triangular factors S.
                              Defined by ( buff_S, ldim_S ).
*/
void __t_geqrf( int t, int m, int n, int skip, fp_t * A, int lda, fp_t * tau, fp_t * T, int ldt, int p) 
{
	A += skip;
	m -= skip;
	int mn = m < n ? m : n;
	
	fp_t *LDWORK = (fp_t*) malloc(sizeof(fp_t) * n * t);

	int k;
	for( k = 0; k < mn; k += t ) {
		int left = mn - k;
		int b = t < left ? t : left;

		int m_ABR = m - k;
		int n_W12 = n - k - b;

		fp_t *buff_ABR = &( A[ lda * k + k ] );          
		fp_t *buff_A12 = &( A[ lda * ( k + b ) + k ] );  
		fp_t *buff_t1  = &( tau[ k ] );                   
		fp_t *buff_T1  = &( T[ t * k + 0 ] );          

		int info;
		LAPACK_geqr2(m_ABR, b, buff_ABR, lda, buff_t1, LDWORK, info );
		if( info != 0 ) {
			printf( "ERROR in NoFLA_QR_tS_unbvar1. Info of dgeqr2: %d\n", info );
		}

		LAPACK_larft(LAPACK_FORWARD, LAPACK_COLWISE, m_ABR, b, buff_ABR, lda, buff_t1, buff_T1, ldt); 

		dlarfb_("Left", "Transpose", "Forward", "Column",\
			&m_ABR, &n_W12, &b,\
			buff_ABR, &lda,\
			buff_T1, &ldt,\
			buff_A12, &lda,\
			LDWORK, &n);
#if 0
		LAPACK_larfb(LAPACK_LEFT, LAPACK_TRANSP, LAPACK_FORWARD, LAPACK_COLWISE,\
			m_ABR, n_W12, b,\
			buff_ABR, lda,\
			buff_T1, ldt,\
			buff_A12, lda,\
			LDWORK, n);
#endif
	}

	free(LDWORK);
}
