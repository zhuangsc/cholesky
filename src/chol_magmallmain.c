#include "chol_magmallmain.h"
#include "chol_kernels_magma.h"
#include "chol_data.h"


void chol_magmall(int mt, int b, int t, REAL **Ah) {
	int k;
	for ( k = 0; k < mt; k++ ) {
		int j;
		for (j=0; j<k; j++) {
			syrk_magmatask(b, Ah[j * mt + k], Ah[k * mt + k], b);
		}

		potrf_magmatask(b, t, Ah[k* mt + k], b );

		int i;
		for (i = k+1; i < mt; i++) {
			int j;
			for (j=0; j<k; j++) {
				gemm_magmatask(b, t, Ah[j * mt + i], Ah[j * mt + k], Ah[ k * mt + i], b);
			}
		}

		for (i = k+1; i < mt; i++) {
	      		trsm_magmatask(b, t, Ah[k * mt + k], Ah[k * mt + i], b);
		}
	}
}
