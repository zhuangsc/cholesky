#ifndef __CHOL_KERNELS_H__
#define __CHOL_KERNELS_H__


#ifdef SINGLE_PRECISION

#define GEMM_TASK 	sgemm_task
#define SYRK_TASK 	ssyrk_task
#define POTRF_TASK 	spotrf_task
#define TRSM_TASK 	strsm_task

#else

#define GEMM_TASK 	dgemm_task
#define SYRK_TASK 	dsyrk_task
#define POTRF_TASK 	dpotrf_task
#define TRSM_TASK 	dtrsm_task

#endif


/* single precision */
#pragma omp task in([b*b]A, [b*b]B) \
		inout ([b*b]C)
void sgemm_task( int b, int t, float *A, float *B, float *C, int ldm);
			      	
#pragma omp task in([b*b]A) \
		inout ([b*b]C) priority(2)
void ssyrk_task(int b, float *A, float *C, int ldm);

#pragma omp task inout([b*b]A) priority(1)
void spotrf_task( int b, int t, float *A, int ldm);

#pragma omp task in([b*b]A) \
		inout ([b*b]B)
void strsm_task( int b, int t, float *A, float *B, int ldm);


/* double precision */
#pragma omp task in([b*b]A, [b*b]B) \
		inout ([b*b]C)
void dgemm_task( int b, int t, double *A, double *B, double *C, int ldm);
			      	
#pragma omp task in([b*b]A) \
		inout ([b*b]C) priority(2)
void dsyrk_task(int b, double *A, double *C, int ldm);

#pragma omp task inout([b*b]A) priority(1)
void dpotrf_task( int b, int t, double *A, int ldm);

#pragma omp task in([b*b]A) \
		inout ([b*b]B)
void dtrsm_task( int b, int t, double *A, double *B, int ldm);


#endif // __CHOL_KERNELS_H__
