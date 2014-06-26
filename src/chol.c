#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>

#include "fptype.h"
#include "chol_config.h"
#include "chol_check.h"
#include "chol_kernels.h"
#include "chol_llmain.h"
#include "chol_rlmain.h"
#include "chol_nllmain.h"
#include "chol_nrlmain.h"
#include "chol_prlmain.h"
#include "chol_setup.h"
#include "chol_utils.h"

#include "chols_llmain.h"


int m; /* order of input matrix */
int ts; /* block size */
int bs; /* sub-block size */
int reps; /* number of repetitions */
int check; /* check result */
int mt; /* number of blocks */
int mr; /* rounded m */
int mtleft;

/* dense */
fp_t *A;
fp_t **Ah;
fp_t *Aorig;	

/* sparse */
void *Ahb;
void *Acsr;
void *Acsc;
int *work;
int format;
long nzA;
long nzL;

//extern pthread_mutex_t mutexhb;
//extern int *vptr_pool, *vpos_pool;
//extern double *vval_pool;
//extern int vptr_unit, vpos_unit, vval_unit;
//extern int vptr_pp, vpos_pp, vval_pp;


int main(int argc, char* argv[]) {
	if ( chol_config(argc, argv) ) {
		return 1;
	}
	
	if ( chol_setup(check, m, mr, ts, bs, mt, mtleft) ) {
		fprintf(stderr, "err: allocating matrix\n");
		return 2;
	}

#if 0
#if USE_PRL 
	octave_matrix_def("chol_octave.m", "A", m, mr, m, mr, mr, ts, Ah[0]);
#else
	octave_matrix_def("chol_octave.m", "A", m, mr, m, mr, ts, ts, Ah[0]);
#endif
#endif

  	unsigned long elapsed= 0 ;
	int r;
	for (r = 0; r < reps; r++) {
		struct timeval start;	
		gettimeofday(&start, NULL);

#if USE_LL
		CHOL_LL(mt, ts, bs, Ah);
#elif USE_RL
		CHOL_RL(mt, ts, bs, Ah);
#elif USE_PRL
		//elapsed += chol_prl(mr, mt, ts, bs, Ah);
#elif USE_NLL
		CHOL_NLL(mt, ts, bs, bs, Ah);
#elif USE_NRL
		CHOL_NRL(mt, ts, bs, bs, Ah);
#elif USE_SLL
		if ( format == MAT_CSC ) {
			chols_ll(Acsc, Acsr, work);
		} else {
			chols_ll_upper(Acsr, Acsc, work);
		}
#endif

		#pragma omp taskwait

		struct timeval stop;	
		gettimeofday(&stop, NULL);

		elapsed += (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec;

		printf("info: numerical fac %lu us\n", elapsed);
  	}

#if 0
#if USE_PRL
	octave_matrix_def("chol_octave.m", "R", m, mr, m, mr, mr, ts, Ah[0]);
#else
	octave_matrix_def("chol_octave.m", "R", m, mr, m, mr, ts, ts, Ah[0]);
#endif
#endif
  

  	int ret = chol_check(check, m, mr, ts, Ah, Aorig);
  

	double elapsedfl = (double) elapsed;
#if 0
  	double dn = (double) m;
  	printf( "time (usec.):\t%.2f GFLOPS: %.4f\n", elapsedfl / (double)reps, (((4.0/3.0)*dn*dn*dn)/((elapsedfl/(double)reps)*1.0e+3)) );
#else
  	FILE *statf=fopen("ompss-stats-0001","w");
  	fprintf(statf,"time : %.2f\n",elapsedfl / (double) reps);
  	fclose(statf);
#endif

#if USE_SLL
	nzL = 0;
	hbmat_t* Z = Acsr;
	for(int i = 0; i < Z->m; ++i){
		for(int j = Z->vptr[i]; j < Z->vptr[i+1]; ++j){
			nzL += ((hbmat_t**)Z->vval)[j]->elemc;
		}
	}

//	printf("nzA %d, nzL %d\n", nzA, nzL);
	if (nzA > nzL)
		nzL = 2 * nzL - Z->m * bs;
	printf("nzA %d nzL %d\n", nzA,nzL);
	statf = fopen("ompss-stats-0001","a");
	fprintf(statf,"degree  : %lf\n", ((double)nzL / (double)nzA));
	fclose(statf);
#endif

	//octave_fin("chol_octave.m", "A", "R");


	chol_shutdown();


	return ret;
}
