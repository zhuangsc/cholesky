#include "chol_setup.h"

#include <stdlib.h>
#include <time.h>

#include "chol_utils.h"
#include "fptype.h"
#include "genmat.h"


extern fp_t *A;
extern fp_t **Ah;
extern fp_t *Aorig;


int chol_setup(int check, int m, int mr, int ts, int bs, int tm, int mleft) { 
  	// Creation of hierarchical matrices. 
	Ah = (fp_t **) malloc( tm * tm * sizeof(fp_t *) );
	if ( Ah == NULL ) {
		return 1;
	}

	// Creation of data matrices
	A = calloc( mr * mr, sizeof(fp_t) );
	if ( A == NULL ) return 1;

	// Initialization of Ah
	int j;
	for ( j=0; j<tm; ++j ) {
    		int i;
		for ( i=0; i<tm; ++i ) {
			Ah[j*tm+i] = (fp_t *) &A[j*ts*mr+i*ts*ts];
		}
  	}

	GENMAT(m, tm, mleft, ts, A);

	Aorig=NULL;
	if ( check ) {
		Aorig = canonical_from_hierarchical( "Aorig", m, mr, m, mr, ts, ts, Ah[0] );
	}

	return 0;
}


void chol_shutdown() {
	free(A);
}
