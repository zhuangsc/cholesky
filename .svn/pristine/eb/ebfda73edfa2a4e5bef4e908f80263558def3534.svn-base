#include "chol_setup.h"
#include "chol_utils.h"
#include <stdlib.h>
#include <time.h>


static double *A;


int chol_setup(int check, int m, int mr, int ts, int bs, int tm, int left, double ***Ah, double **Aorig) {
	// Creation of hierarchical matrices. 
	*Ah = (double **) malloc( tm * sizeof(double *) );
	if(*Ah == NULL) return 1;

	double **lAh = *Ah;

	// Creation of data matrices
	A = (double *) calloc( mr * mr, sizeof(double) );
	if (A==NULL) return 1;

	// Initialization of Ah
	int j;
	for (j = 0; j < tm; j++) {
		lAh[j] = (double *) &A[ts*mr * j];
  	}

	// Initialization of A with random values
	srand48( ( long int ) time(NULL) );

	double *Ain = &A[0];
	for (j = 0; j < mr; j++) {
		int i;
		for( i = 0; i < mr; i++) {
			if(i<m && j<m) {
				double dran = drand48();
				if(dran == 0.0)  printf("generated 0\n");
				if(j==i) {
					dran += m;
				}
				*Ain = dran;
			} else {
				if(i==j) {
					*Ain=1.00;
				}
			}
			//printf("(%i,%i) = %.2f\n",i,j,*Ain);
			Ain++;
		}
  	}

	*Aorig=NULL;
	if(check) {
		*Aorig = canonical_from_hierarchical( "Aorig", m, mr, m, mr, mr, ts, lAh[0] );
	}

	return 0;
}


void chol_shutdown() {
	free(A);
}
