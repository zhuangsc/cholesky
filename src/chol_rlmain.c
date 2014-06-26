#include "chol_rlmain.h"

#include <sys/time.h>
#include <stddef.h>

#include "fptype.h"
#include "chol_kernels.h"


int CHOL_RL(int mt, int b, int t, fp_t **Ah) {
	int k;
	for ( k = 0; k < mt ;k++ ) {

		POTRF_TASK(b, t, Ah[k* mt + k], b);

		int i;
		for (i = k+1; i < mt; i++) {
		      	TRSM_TASK(b, t, Ah[k * mt + k], Ah[k * mt + i], b);
		}

		for (i = k+1; i < mt; i++) {
			int j;
			for (j=k+1; j<i; j++) {
				GEMM_TASK(b, t, Ah[k * mt + i], Ah[k * mt + j], Ah[j * mt + i], b);
			}
			SYRK_TASK(b, Ah[k * mt + i], Ah[i * mt + i], b);
		}
	}

	return 0;
}
