#ifndef __BLAS_EXT_H__
#define __BLAS_EXT_H__


#if SINGLE_PRECISION
#define BLAS_set 	BLAS_sset
#else
#define BLAS_set 	BLAS_dset
#endif


static inline __attribute__((always_inline)) void BLAS_dset(double val, int m_A, int n_A, double * A, int ldim_A)
{
  int  i, j;
  for( j = 0; j < n_A; j++ ) {
    for( i = 0; i < m_A; i++ ) {
      A[ ldim_A * j + i ] = val;
    }
  }
}


static inline __attribute__((always_inline)) void BLAS_sset(float val, int m_A, int n_A, float * A, int ldim_A)
{
  int  i, j;
  for( j = 0; j < n_A; j++ ) {
    for( i = 0; i < m_A; i++ ) {
      A[ ldim_A * j + i ] = val;
    }
  }
}


// B := A 
static inline __attribute__((always_inline)) int NoFLA_Copy( int m_A, int n_A, double * A, int ldim_A, double * B, int ldim_B ) 
{
  int  i, j;

  for( j = 0; j < n_A; j++ ) {
    for( i = 0; i < m_A; i++ ) {
      B[ ldim_B * j + i ] = A[ ldim_A * j + i ];
    }
  }
  return 0;
}


// triu( B ) := triu( A )
static inline __attribute__((always_inline)) int NoFLA_Copy_triu( int m_A, int n_A, double * A, int ldim_A, double * B, int ldim_B ) 
{
  int     i, j;

  for ( j = 0; j < n_A; j++ ) {
    for ( i = 0; i <= min( j, m_A - 1 ); i++ ) {
      B[ ldim_B * j + i ] = A[ ldim_A * j + i ];
    }
  }
  return 0;
}


// B := alpha * A + B 
static inline __attribute__((always_inline)) int NoFLA_Axpy( int m_A, int n_A, double alpha, double * A, int ldim_A, double * B, int ldim_B ) 
{
  int  i, j;
  for( j = 0; j < n_A; j++ ) {
    for( i = 0; i < m_A; i++ ) {
      B[ ldim_B * j + i ] += alpha * A[ ldim_A * j + i ];
    }
  }
  return 0;
}


// Zero strict lower part 
static inline __attribute__((always_inline)) int NoFLA_Zero_strict_lower_part( int m_A, int n_A, double * A, int ldim_A ) 
{
  int  i, j;

  for( j = 0; j < n_A; j++ ) {
    for( i = j+1; i < m_A; i++ ) {
      A[ ldim_A * j + i ] = 0.0;
    }
  }
  return 0;
}


#endif // __BLAS_EXT_H__
