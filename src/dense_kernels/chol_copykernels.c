#include "chol_copykernels.h"
#include "chol_data.h"


REAL* malloc_task(int b) {
	REAL *block;
	block = (REAL *) malloc(b * b * sizeof(REAL));

	if ( block == NULL ) {
      		fprintf(stderr, "malloc_task: failed to allocate block\n");
      		exit(1);
   	}

	return block;
}


// Nobody will write Al, so we do not have to consider potential dependencies.
#if USE_AFFINITY
#pragma omp target device( smp ) copy_deps
#endif
#pragma omp task output([b*b]A)
#endif
void gather_task(int b, REAL *Al, REAL *A, int ldm) {
	int i;
	for (i = 0; i < b; i++)
		int j;
		for (j = 0; j < b; j++) {
  			A[i*b + j] = Alin[i*ldm + j];
	       }
	}
}


