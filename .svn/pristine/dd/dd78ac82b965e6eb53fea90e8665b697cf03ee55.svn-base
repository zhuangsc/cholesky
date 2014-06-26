#include "chols_kernels.h"

#include <math.h>

#include "array.h"


/* A in CSR, lower-triangular */
void potrf_sparse_upper(hbmat_t* A) {
 	/* Check if the input matrix has properly set */
	if ( A->vptr == NULL ) {
		hyper_sym_csr_task2(A);
	}

	int m = A->m;
	int* vptr = A->vptr; int* vpos = A->vpos; 
	double* vval = A->vval;

	int d;
	for ( d=0 ; d<m; ++d ) {
		int dstart = vptr[d];
		int dend = vptr[d+1];

		int diagidx = dend - 1;
		double diagel = vval[diagidx];
		/* compute element on diagonal (d,d) */
		int j;
		for ( j = dstart; j < diagidx; ++j ) {
			diagel -= vval[j]*vval[j];
		}
		diagel = sqrt(diagel);
		vval[diagidx] = diagel;


		int i;
		for ( i = d+1; i < m; ++i ) {
			int jx = vptr[i];
			int j = vpos[jx];
			int jend = vptr[i+1];

			int djx = dstart;
			int dj = vpos[djx];

			double acc = 0;
			while ( dj < d && djx < dend && j < d && jx < jend ) {
				int djadv = dj < j; 
				int jadv  = j  < dj; 

				jx  += jadv? 1 : 0;
				djx += djadv? 1 : 0;
				dj = vpos[djx];
			 	j  = vpos[ jx];

				if ( dj == j && j < d && djx < dend && jx < jend ) {
					acc += vval[djx++] * vval[jx++];
					dj = vpos[djx];
			 		j  = vpos[ jx];
				}
			}

			while ( vpos[jx] < d && jx < jend ) { ++jx; }

			if ( jx < jend && vpos[jx]==d ) {
				double subdval = vval[jx] - acc;
				vval[jx] = subdval / diagel;
			}
		}
	}
}


/* A and C in CSR, lower-triangular */
void dsyrk_sparse_upper(hbmat_t* A, hbmat_t* C) {
	 /* Check if the input matrix has properly set */
	if ( A->vptr == NULL )
		hyper_sym_csr_task2(A);
	if ( C->vptr == NULL )
		hyper_sym_csr_task2(C);

	int n = A->n; int m = A->m;
	int* vptr = A->vptr; int* vpos = A->vpos; double* vval = A->vval;
	int* vptr_c = C->vptr; int* vpos_c = C->vpos; double* vval_c = C->vval;

	double* peela = malloc(2 * m * sizeof(double));
	double* peelc = &peela[m];

	int i;
	for ( i = 0; i < m; ++i ) {
//		array_clear(peelc, m);
		array_s2d(A, peela, i);

		mkl_cspblas_dcsrgemv("N", &m, vval, vptr, vpos, peela, peelc);

		int k;
		for ( k = vptr_c[i]; k < vptr_c[i+1]; ++k ) {
			int col_pos = vpos_c[k];
			vval_c[k] -= peelc[col_pos];
		}

		array_clear2(A, peela, i);
	}

	free(peela); 
	//free(peelc);
}


/* A, B and C in CSR */
/* C = C - A * B^T, or */
/* C^T = C^T - B * A^T */
void dgemm_sparse_upper(hbmat_t* A, hbmat_t* B, hbmat_t* C) {
 	/* Check if the input matrix has properly set */
	if ( A->vptr == NULL ) {
		hyper_sym_csr_task2(A);
	}

//	if ( B->vptr == NULL )
//		hyper_sym_csr_task1(B);

	if ( C->vptr == NULL ) {
		hyper_sym_csr_task2(C);
	}

	int m = B->m; int n = B->n;
	int* vptr = B->vptr; int* vpos = B->vpos; double* vval = B->vval;
	int* vptr_c = C->vptr; int* vpos_c = C->vpos; double* vval_c = C->vval;

	double* peelb = malloc(2*m*sizeof(double));
	double* peelc = &peelb[m];

	int i;
	for ( i = 0; i < m; ++i ) {
		array_clear(peelb, m);
		array_s2d(A, peelb, i);
		mkl_cspblas_dcsrgemv("N", &m, vval, vptr, vpos, peelb, peelc);
		int k;
		for ( k = vptr_c[i]; k < vptr_c[i+1]; k++ ) {
			int col_pos = vpos_c[k];
			vval_c[k] -= peelc[col_pos];
		}
	}

	free(peelb); 
}


/* A CSC, C CSC */
/* C = C \ A */
void dtrsm_sparse_upper(hbmat_t* A, hbmat_t* C){
	/* Check if the input matrix has properly set */
#if 0
	if ( A->vptr == NULL )
		hyper_sym_csr_task1(A);
#endif

	if ( C->vptr == NULL ) {
		hyper_sym_csr_task2(C);
	}

	int n = A->n; int m = A->m;
	int* vptr = A->vptr; int* vpos = A->vpos; double* vval = A->vval;

	double* peelb = calloc(2*m, sizeof(double));
	double* peelc = &peelb[m];
	double alpha = 1;

	int i;
	for ( i = 0; i < n; ++i ) {
		array_s2d(C, peelb, i);
		mkl_dcsrsv("N", &n, &alpha, "TLNC", vval, vpos, vptr, vptr+1, peelb, peelc);
		array_d2s(C, peelc, i);
		//array_clear(peelb, m);
		array_clear2(C, peelb, i);
	}

	free(peelb); 
}
