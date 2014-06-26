#if OMPSS_HYBRID || OMPSS_CUDA


#include "chol_kernels.h"
#include "chol_kernels_cuda.h"
#include "chol_data.h"
#include "bsc_potrf.h"

#include <cublas.h>

void gemm_cudatask( int b, int t, REAL *A, REAL *B, REAL *C, int ldm) {
	REAL mone = -1.0; 
	REAL one = 1.0;

	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(nanos_get_cublas_handle(), stream);
	cublasHandle_t handle = nanos_get_cublas_handle();

	CUBLASGEMM(handle, CUBLAS_OP_N, CUBLAS_OP_T, b, b, b, &mone, A, ldm,
		B, ldm, &one, C, ldm);
}
			      	

void syrk_cudatask(int b, REAL *A, REAL *C, int ldm) {
	REAL mone = -1.0;
	REAL one = 1.0;

	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(nanos_get_cublas_handle(), stream);
	cublasHandle_t handle = nanos_get_cublas_handle();

	CUBLASSYRK(handle, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, b, b, &mone, A, ldm, &one, C, ldm);
}


void potrf_cudatask(int b, int t, REAL *A, int ldm) {
	cublasHandle_t handle = nanos_get_cublas_handle();
	int info;
	bsc_potrf( handle, 'L', b, A, ldm, &info );
}


void trsm_cudatask( int b, int t, REAL *A, REAL *B, int ldm) {
    	REAL one = 1.0;

	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetStream(nanos_get_cublas_handle(), stream);
	cublasHandle_t handle = nanos_get_cublas_handle();

	CUBLASTRSM(handle, CUBLAS_SIDE_RIGHT, CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, b, b, &one, A, ldm, B, ldm );
}


#undef CUBLASGEMM
#undef CUBLASSYRK
#undef CUBLASTRSM


#endif // OMPSS_HYBRID || OMPSS_CUDA
