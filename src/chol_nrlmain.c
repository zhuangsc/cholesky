#include "chol_nrlmain.h"

#include "fptype.h"
#include "chol_kernels_nested.h"


int CHOL_NRL(int mt, int sb, int b, int t, fp_t **Ah) {
	//printf("chol_nrl: %i %i %i %i\n", mt, sb, b, t);
	int k;
	for ( k = 0; k < mt ;k++ ) {

		POTRF_NTASK(sb, b, t, Ah[k* mt + k]);

		int i;
		for (i = k+1; i < mt; i++) {
		      	TRSM_NTASK(sb, b, t, Ah[k * mt + k], Ah[k * mt + i]);
		}

		for (i = k+1; i < mt; i++) {
			int j;
			for (j=k+1; j<i; j++) {
				GEMM_NTASK(sb, b, t, Ah[k * mt + i], Ah[k * mt + j], Ah[j * mt + i]);
			}
			SYRK_NTASK(sb, b, t, Ah[k * mt + i], Ah[i * mt + i]);
		}
	}

	return 0;
}
