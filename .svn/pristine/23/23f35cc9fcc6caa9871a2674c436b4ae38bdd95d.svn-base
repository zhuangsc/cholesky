#include "ompss_cholesky.h"

#include "chol_llmain.h"
#include "chol_rlmain.h"
#include "chol_nllmain.h"
#include "chol_nrlmain.h"
#include "chols_llmain.h"

#include "hb.h"
#include "hbext.h"
#include "hbconvrt.h"

#include <stdio.h>


/* hypermatrix */
int ompss_dchol_hll(int mt, int b, int t, double **Ah) {
	return dchol_ll(mt, b, t, Ah);
}

int ompss_schol_hll(int mt, int b, int t, float **Ah) {
	return schol_ll(mt, b, t, Ah);
}

int ompss_dchol_hrl(int mt, int b, int t, double **Ah) {
	return dchol_rl(mt, b, t, Ah);
}

int ompss_schol_hrl(int mt, int b, int t, float **Ah) {
	return schol_rl(mt, b, t, Ah);
}


/* unimplemented */
int ompss_dchol_ll(int n, int b, int t, double *A) {
	fprintf(stderr, "err: not yet implemented\n");
	return 0;
}

int ompss_schol_ll(int n, int b, int t, float *A) {
	fprintf(stderr, "err: not yet implemented\n");
	return 0;
}

int ompss_dchol_rl(int n, int b, int t, double *A) {
	fprintf(stderr, "err: not yet implemented\n");
	return 0;
}

int ompss_schol_rl(int n, int b, int t, float *A) {
	fprintf(stderr, "err: not yet implemented\n");
	return 0;
}


/* nested */
int ompss_dchol_nhll(int mt, int b, int t, double **Ah) {
	return dchol_nll(mt, b, t, t, Ah);
}

int ompss_schol_nhll(int mt, int b, int t, float **Ah) {
	return schol_nll(mt, b, t, t, Ah);
}

int ompss_dchol_nhrl(int mt, int b, int t, double **Ah) {
	return dchol_nrl(mt, b, t, t, Ah);
}

int ompss_schol_nhrl(int mt, int b, int t, float **Ah) {
	return schol_nrl(mt, b, t, t, Ah);
}



/* sparse, 0-based hbmat_t */
hbmat_t* ompss_csr_dchol_ll(int b, hbmat_t *A, int *work) {

	int *tree = etree(A);
	hbmat_t *Acsr = hb2hbh_hyper_sym_csr(A, b, tree);

	hbmat_t *Acsc = Acsr->trans = hb2csr(Acsr);

	chols_ll_upper(Acsr, Acsc, work);

	return Acsr;
}

hbmat_t* ompss_csc_dchol_ll(int b, hbmat_t *A, int *work) {
	hbmat_t *tAcsr = A->trans = hb2csr(A);
	int *tree = etree(tAcsr);
	free(tAcsr);

	hbmat_t *Acsc = hb2hbh_sym_etree(A, b, tree);
	hbmat_t *Acsr = Acsc->trans = hb2csr(Acsc);

	free(tree);

	chols_ll(Acsc, Acsr, work);

	return Acsc;
}

