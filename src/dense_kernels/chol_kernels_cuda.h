#if defined(OMPSS_HYBRID) || defined(OMPSS_CUDA)


#ifndef __CHOL_KERNELS_CUDA_H__
#define __CHOL_KERNELS_CUDA_H__


#include "chol_kernels.h"
#include "chol_data.h"


#pragma omp target device (cuda) implements (gemm_task) copy_deps
#pragma omp task in ([b*b]A, [b*b]B) \
		inout ([b*b]C)
void gemm_cudatask( int b, int t, REAL *A, REAL *B, REAL *C, int ldm);
			      	

#pragma omp target device (cuda) implements (syrk_task) copy_deps
#pragma omp task in ([b*b]A) \
		inout ([b*b]C) priority(1)
void syrk_cudatask(int b, REAL *A, REAL *C, int ldm);


#pragma omp target device (cuda) implements (potrf_task) copy_deps
#pragma omp task inout([b*b]A) priority(3)
void potrf_cudatask(int b, int t, REAL *A, int ldm);


#pragma omp target device (cuda) implements (trsm_task) copy_deps
#pragma omp task in([b*b]A) \
		inout ([b*b]B) priority(2)
void trsm_cudatask( int b, int t, REAL *A, REAL *B, int ldm);


#endif // __CHOL_KERNELS_CUDA_H__


#endif // OMPSS_HYBRID || OMPSS_CUDA
