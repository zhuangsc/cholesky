#include "chol_nllmain.h"

#include "fptype.h"
#include "chol_kernels_nested.h"


int CHOL_NLL(int mt, int sb, int b, int t, fp_t **Ah) {
	int k;
	for ( k = 0; k < mt; k++ ) {
		int j;
		for ( j=0; j<k; j++ ) {
		//	#pragma omp task in([sb*sb]Ah) inout([sb*sb]Ah) priority( (mt-j)+10 ) untied
			SYRK_NTASK(sb, b, t, Ah[j * mt + k], Ah[k * mt + k]);
		}

 		//#pragma omp task inout([sb*sb]Ah) priority( 100000 ) untied
		POTRF_NTASK(sb, b, t, Ah[k* mt + k]);

		int i;
		for ( i = k+1; i < mt; i++ ) {
			int j;
			for ( j=0; j<k; j++ ) {
				GEMM_NTASK(sb, b, t, Ah[j * mt + i], Ah[j * mt + k], Ah[k * mt + i]);
			}

			TRSM_NTASK(sb, b, t, Ah[k * mt + k], Ah[k * mt + i]);
		}

#if 0
		for (i = k+1; i < mt; i++) {
			#pragma omp task in([sb*sb]Ah[k*mt + k]) inout([sb*sb]Ah[k*mt + i]) priority( (mt-i)+10 ) untied
			trsm_ntask(sb, b, t, Ah[k * mt + k], Ah[k * mt + i]);
		}
#endif
	}

	return 0;
}
