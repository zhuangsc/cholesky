#include "chols_llmain.h"

#include "chols_kernels.h"
#include "hb.h"
#include "hbconvrt.h"
#include "matfprint.h"

#include <sys/time.h>

/*
 * CSC
 */
int chols_ll(hbmat_t* A_csc, hbmat_t* A_csr, int *work) {
	int n = A_csc->n;

	int* vptr_csc = A_csc->vptr; int* vpos_csc = A_csc->vpos;
	hbmat_t** vval_csc = A_csc->vval;

	int* vptr_csr = A_csr->vptr; int* vpos_csr = A_csr->vpos;
	hbmat_t** vval_csr = A_csr->vval;

	for ( int J=0; J<n; ++J ) {
		hbmat_t* l22 = vval_csc[vptr_csc[J]];
		for ( int k = vptr_csr[J]; k < vptr_csr[J+1]-1; ++k ) {
			dsyrk_sparse(vval_csr[k], l22);
		}

		potrf_sparse(l22);

		int acc = 0;
		for ( int b = vptr_csc[J]+1; b < vptr_csc[J+1]; ++b, ++acc ) {
			work[acc] = vptr_csr[vpos_csc[b]];
		}

		for ( int j = vptr_csr[J]; j < vptr_csr[J+1]-1; ++j ) {
			int col_21 = vpos_csr[j];
			acc = 0;
			for ( int k = vptr_csc[J]+1; k < vptr_csc[J+1]; ++k, ++acc ) {
				int row_32 = vpos_csc[k];
				int cacc = work[acc];
				int col_32 = vpos_csr[cacc];
				while ( col_32 < col_21 && cacc < vptr_csr[row_32+1]) {
					col_32 = vpos_csr[++cacc];
				}
				if ( col_32 == col_21 ) {
					dgemm_sparse(vval_csr[cacc], vval_csr[j], vval_csc[k]);
				}
				work[acc] = cacc;
			}
		}

		for ( int j = vptr_csc[J]+1; j < vptr_csc[J+1]; ++j ) {
			dtrsm_sparse(l22, vval_csc[j]);
		}
	}

#if 1
#pragma omp taskwait
	hb_free(A_csr);
	hbmat_t* test = hbh2hb(A_csc);
	hb_print_CSC2("L000.dat", test);
	hb_free(test);
#endif

	return 0;
}

/* A_csr: CSR, lower-triangular, A_csc: CSC, lower-trangular */
int chols_ll_upper(hbmat_t* A_csr, hbmat_t* A_csc, int *work) {
	int n = A_csr->n;

	int* vptr_csc = A_csc->vptr; int* vpos_csc = A_csc->vpos;
	hbmat_t** vval_csc = A_csc->vval;

	int* vptr_csr = A_csr->vptr; int* vpos_csr = A_csr->vpos;
	hbmat_t** vval_csr = A_csr->vval;

	for ( int d=0; d<n; ++d ) {
		/* diagonal block D */
		hbmat_t* D = vval_csc[vptr_csc[d]];

		/* symmetric rank-k update of D with the blocks on the left */
		for ( int k = vptr_csr[d]; k < vptr_csr[d+1]-1; ++k ) {
			dsyrk_sparse_upper(vval_csr[k], D);
		}

		/* factorize D */
		potrf_sparse_upper(D);

		int acc = 0;
		for ( int b = vptr_csc[d]+1; b < vptr_csc[d+1]; ++b, ++acc ) {
			work[acc] = vptr_csr[vpos_csc[b]];
		}

		/* update the blocks underneath D */
		for ( int j = vptr_csr[d]; j < vptr_csr[d+1]-1; ++j ) {
			int col_21 = vpos_csr[j];
			int acc = 0;

			for ( int k = vptr_csc[d]+1; k < vptr_csc[d+1]; ++k, ++acc ) {
				int row_32 = vpos_csc[k];
				int cacc = work[acc];
				int col_32 = vpos_csr[cacc];

				while ( col_32 < col_21 && cacc < vptr_csr[row_32+1]) {
					col_32 = vpos_csr[++cacc];
				}
				if ( col_32 == col_21 ) {
					dgemm_sparse_upper(vval_csr[cacc], vval_csr[j], vval_csc[k]);
				}
				work[acc] = cacc;
			}
		}

		for ( int j = vptr_csc[d]+1; j < vptr_csc[d+1]; j++ ) {
			dtrsm_sparse_upper(D, vval_csc[j]);
		}
	}


#if 0
#pragma omp taskwait
	hb_free(A_csc);
	hbmat_t* test = hbh2hb(A_csr);
	fprint_csc("L000.dat", test, OFFS_C);
	hb_free(test);
#endif

	return 0;
}
