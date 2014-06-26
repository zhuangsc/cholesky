#ifndef __CHOL_KERNELS_NESTED_H__
#define __CHOL_KERNELS_NESTED_H__


#ifdef SINGLE_PRECISION

#define GEMM_NTASK 		sgemm_ntask
#define SYRK_NTASK 		ssyrk_ntask
#define POTRF_NTASK 	spotrf_ntask
#define TRSM_NTASK 		strsm_ntask

#else

#define GEMM_NTASK 		dgemm_ntask
#define SYRK_NTASK 		dsyrk_ntask
#define POTRF_NTASK 	dpotrf_ntask
#define TRSM_NTASK 		dtrsm_ntask

#endif


/* single precision */
#pragma omp target device (smp) 
#pragma omp task in ([m*m]A, [m*m]B) \
		inout ([m*m]C)
void sgemm_ntask( int m, int b, int t, float *A, float *B, float *C); 
			      	

#pragma omp target device (smp) 
#pragma omp task in ([m*m]A) \
		inout ([m*m]C) priority(2)
void ssyrk_ntask(int m, int b, int t, float *A, float *C);


#pragma omp target device (smp) 
#pragma omp task inout([m*m]A) priority(1)
void spotrf_ntask( int m, int b, int t, float *A );


#pragma omp target device (smp) 
#pragma omp task in([m*m]A) \
		inout ([m*m]B)
void strsm_ntask( int m, int b, int t, float *A, float *B);


/* single precision */
#pragma omp target device (smp) 
#pragma omp task in ([m*m]A, [m*m]B) \
		inout ([m*m]C)
void dgemm_ntask( int m, int b, int t, double *A, double *B, double *C); 
			      	

#pragma omp target device (smp) 
#pragma omp task in ([m*m]A) \
		inout ([m*m]C) priority(2)
void dsyrk_ntask(int m, int b, int t, double *A, double *C);


#pragma omp target device (smp) 
#pragma omp task inout([m*m]A) priority(1)
void dpotrf_ntask( int m, int b, int t, double *A );


#pragma omp target device (smp) 
#pragma omp task in([m*m]A) \
		inout ([m*m]B)
void dtrsm_ntask( int m, int b, int t, double *A, double *B);


#endif // __CHOL_KERNELS_NESTED_H__
