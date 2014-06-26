#ifndef __QRCA_UTILS_H__
#define __QRCA_UTILS_H__


#include <stdio.h>


#define min( a, b ) ( (a) < (b) ? (a) : (b) )
#define max(a,ts) ((a) >= (ts) ? (a) : (ts))


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



void print_matrix_blocked( FILE* str, char * name, int m, int mr, int n, int nr, int br, int bc, double * A );
double * canonical_from_hierarchical ( char * name, int m, int mr, int n, int nr, int br, int bc, double * A );
void print_matrix_canonical( char *name, int m, int n, double *A, int lda );


#endif // __QRCA_UTILS_H__
