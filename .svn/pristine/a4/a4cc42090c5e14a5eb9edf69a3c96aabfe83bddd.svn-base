#include "chol_prlmain.h"
#include "chol_pkernels.h"

#include <sys/time.h>
#include <stddef.h>


unsigned long chol_prl(int m, int mt, int b, int t, double **Ah) {
	struct timeval start, stop;
	
	gettimeofday(&start, NULL);

	int skipa = 0;
	//printf("m %i mt %i b %i t %i\n", m, mt, b, t);

	int k;
	for ( k = 0; k < mt; k++ ) {
		fac_task(skipa, m, b, t, Ah[ k ]);

		int skipb = skipa;
		int j;
		for (j=k+1; j<mt; j++) {
			skipb += b;
			upd_task(skipa, skipb, m, b, t, Ah[k], Ah[j]);
		}
		skipa += b;
	}

#pragma omp taskwait

	gettimeofday(&start, NULL);

	return (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec;
}
