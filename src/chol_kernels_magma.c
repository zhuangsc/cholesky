#include "chol_kernels.h"
#include "chol_data.h"


void gemm_magmatask( int m, int t, REAL *A, REAL *B, REAL *C) {
	REAL mone = -1.0; 
	REAL one = 1.0;

	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetKernelStream(stream);

 	MAGMABLASGEMM('N', 'T', m, m, m, mone, A, m, B, m, one, C, m);
}
			      	

void syrk_magmatask(int m, REAL *A, REAL *C) {
	REAL mone = -1.0;
	REAL one = 1.0;

	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetKernelStream(stream);

 	MAGMABLASSYRK('L', 'N', m, m, mone, A, m, one, C, m );
}


void potrf_magmatask( int m, int t, REAL *A ) {
 	//cudaStream_t stream = nanos_get_kernel_execution_stream();
	//cublasSetKernelStream(stream);

	int info;
	MAGMAPOTRF('L', m, A, m, &info );
}


void trsm_magmatask( int m, int t, REAL *A, REAL *B) {
    	REAL one = 1.0;
 	cudaStream_t stream = nanos_get_kernel_execution_stream();
	cublasSetKernelStream(stream);

	MAGMABLASTRSM('R', 'L', 'T', 'N', m, m, one, A, m, B, m );
}


#undef MAGMABLASGEMM
#undef MAGMABLASSYRK
#undef MAGMAPOTRF
#undef MAGMABLASTRSM
