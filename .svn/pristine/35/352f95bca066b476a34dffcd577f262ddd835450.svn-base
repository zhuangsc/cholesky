#include "chol_kernels_nested.h"

#include "fptype.h"
#include "chol_kernels.h"



void GEMM_NTASK( int m, int b, int t, fp_t *A, fp_t *B, fp_t *C) {
	int k;
	for ( k=0; k<m ;k+=b) {
		int i;
		for(i=0; i<m; i+=b) {
			int j;
			for (j=0; j<m; j+=b) {
				GEMM_TASK(b, t, &A[k*m+i], &B[k*m+j], &C[j*m+i], m);    
			}
		}
   	}

	#pragma omp taskwait
}
			      	

void SYRK_NTASK(int m, int b, int t, fp_t *A, fp_t *C) {
	int k;
	for (k=0; k<m; k+=b) {
		int j;
		for (j=0; j<m; j+=b) {
			SYRK_TASK(b, &A[j*m + k], &C[k*m + k], m);   
        	}
        
		int i;
		for (i=0; i<k; i+=b) {
			int j;
			for (j=0; j<m; j+=b) {
               			GEMM_TASK(b, t, &A[j*m+k], &A[j*m+i], &C[i*m+k], m); 
            		}
        	}
	}
        
	#pragma omp taskwait
}


void POTRF_NTASK( int m, int b, int t, fp_t *A ) {
#if 0
	nanos_wd_t wd = nanos_current_wd();
	unsigned priority = nanos_get_wd_priority( wd );
	if ( priority != 100000 ) {
       		fprintf( stderr, "nested potrf: priority != 100000\n" );
	}
#endif

	int k;
	for (k = 0; k < m; k += b) {
		POTRF_TASK(b, t, &A[k*m+k], m);               

		int i;
		for (i = k + b; i < m; i+=b) {                     
			TRSM_TASK(b, t, &A[k*m + k], &A[k*m + i], m);
       		}

		for (i = k + b; i < m; i+=b) {                   
			int j;
  			for (j = k + b; j < i; j+=b) {
				GEMM_TASK(b, t, &A[k*m + i], &A[k*m + j], &A[j*m + i], m);
          		}
	  		SYRK_TASK(b, &A[k*m + i], &A[i*m + i], m);
       		}
    	}

	#pragma omp taskwait
}


void TRSM_NTASK( int m, int b, int t, fp_t *A, fp_t *B) {
	int k;
	for (k=0; k<m; k+=b) {
		int i;
		for (i=0; i<m ; i+=b) {
			TRSM_TASK(b, t, &A[k*m + k], &B[k*m + i], m);
			int j;
			for (j=k+b; j<m; j+=b) {
				GEMM_TASK(b, t, &B[k*m + i], &A[k*m + j], &B[j*m + i], m);
			}
		}
	}

	#pragma omp taskwait
}
