#include "hb.h"


/* Assumes A is square, pad for Cholesky factorization */
hbmat_t * csc_pad(hbmat_t *A, int mr, int offs) {
	int m = A->m;
	if ( mr <= m ) {
		return A;
	}

	int elemc = A->elemc;

	int *vptr = A->vptr; int *vpos = A->vpos;
	double *vval = (double*) A->vval;

	int mdif = mr - m;
	int mrp1 = mr + 1;
	int nelemc = elemc + mdif;

	vptr = realloc(vptr, sizeof(int) * mrp1 );
	vpos = realloc(vpos, sizeof(int) * nelemc );
	vval = realloc(vval, sizeof(double) * nelemc );
	if ( vptr == NULL || vpos == NULL ||  vval == NULL ) {
		fprintf(stderr, "pad: could not allocate for padding\n");
		return NULL;
	}

	int mp1 = m + 1;
	int i;
	for ( i=mp1; i<mrp1; ++i ) {
		vptr[i] = vptr[i-1] + 1;
	}
	for ( i=elemc; i<nelemc; ++i ) {
		vpos[i] = m++ + offs;
		vval[i] = 1.0;
	}

	A->vptr = vptr;
	A->vpos = vpos;
	A->vval = vval;
	A->elemc = nelemc;
	A->m = mr;
	A->n = mr;

	return A;
}
