#ifndef __CHOLS_KERNELS_H__
#define __CHOLS_KERNELS_H__


#include "hb.h"


#pragma omp task inout([1]A) priority(1)
void potrf_sparse(hbmat_t* A);

#pragma omp task in([1]A) inout([1]B) priority(2)
void dsyrk_sparse(hbmat_t* A, hbmat_t* B);

#pragma omp task in([1]A, [1]B) inout([1]C)
void dgemm_sparse(hbmat_t* A, hbmat_t* B, hbmat_t* C);

#pragma omp task in([1]A) inout([1]B)
void dtrsm_sparse(hbmat_t* A, hbmat_t* B);

#pragma omp task inout([1]A) priority(4)
void potrf_sparse_upper(hbmat_t* A);

#pragma omp task in([1]A) inout([1]B) priority(3)
void dsyrk_sparse_upper(hbmat_t* A, hbmat_t* B);

#pragma omp task in([1]A, [1]B) inout([1]C)
void dgemm_sparse_upper(hbmat_t* A, hbmat_t* B, hbmat_t* C);

#pragma omp task in([1]A) inout([1]B)
void dtrsm_sparse_upper(hbmat_t* A, hbmat_t* B);


#endif
