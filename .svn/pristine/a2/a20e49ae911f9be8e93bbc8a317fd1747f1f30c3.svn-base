#include "chol_llmain.h"

#include <sys/time.h>
#include <stddef.h>

#include "chol_kernels.h"
#include "fptype.h"


int CHOL_LL(int mt, int b, int t, fp_t **Ah) {
	int k;
	for ( k = 0; k < mt; k++ ) {
		int j;
		for (j=0; j<k; j++) {
			SYRK_TASK(b, Ah[j * mt + k], Ah[k * mt + k], b);
		}

		POTRF_TASK(b, t, Ah[k* mt + k], b);

		int i;
		for (i = k+1; i < mt; i++) {
			int j;
			for (j=0; j<k; j++) {
				GEMM_TASK(b, t, Ah[j * mt + i], Ah[j * mt + k], Ah[ k * mt + i], b);
			}

		  	TRSM_TASK(b, t, Ah[k * mt + k], Ah[k * mt + i], b);
		}

#if 0
		for (i = k+1; i < mt; i++) {
		      TRSM_TASK(b, t, Ah[k * mt + k], Ah[k * mt + i] );
		}
#endif
	}

	return 0;
}
