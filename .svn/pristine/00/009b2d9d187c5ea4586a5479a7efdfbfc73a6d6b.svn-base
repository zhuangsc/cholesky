#include "chols_kernels.h"


void potrf_sparse_upper(hbmat_t* A) {

	/*
	 * Check if the input matrix has properly set
	 */
	if ( A->vptr == NULL )
		hyper_sym_csr_task2(A);

	int m = A->m;
	int* vptr = A->vptr; int* vpos = A->vpos; 
	double* vval = A->vval;

	int I;
	for ( I = 0 ; I < m; I++) {
		//l22 = sqrt (a22 - l12*l12_t)
		double* l22 = &(vval[vptr[I+1]-1]);
		int p_a22 = vpos[vptr[I+1]-1];

		int j;
		for ( j = vptr[I]; j < vptr[I+1]-1; j++ ) {
			*l22 -= vval[j]*vval[j];
		}
		*l22 = sqrt(*l22);
		
		// l32 = (a32-L32*l12)/l22
		double l31l12 = 0;
		int i;
		for ( i = I+1; i < m; i++){
			int j;
			for( j = vptr[i]; j < vptr[i+1]; j++) {
				if (vpos[j] < I)
					continue;

				if (vpos[j] == I){
					double *l32 = &(vval[j]);
					l31l12=0;
					int p_a32 = vptr[i];
					int p_a33 = vptr[i+1];
					int k;
					for( k = vptr[I]; k < vptr[I+1]-1; k++) {
						int l;
						for ( l = p_a32; l < p_a33; l++) {
							if( vpos[l] < vpos[k] ){
								p_a32 = l;
								continue;
							}
							if (vpos[l] == vpos[k]){
								l31l12 += vval[l] * vval[k];
								break;
							}
							if (vpos[l] > vpos[k])
								break;
						}
					}
					*l32 = (*l32 - l31l12)/(*l22);
					break;
				}

				if (vpos[j] > I)
					break;
			}
		}
	}
}

void dsyrk_sparse_upper(hbmat_t* A, hbmat_t* C){

	/*
	 * Check if the input matrix has properly set
	 */
	if ( A->vptr == NULL )
		hyper_sym_csr_task2(A);

	if ( C->vptr == NULL )
		hyper_sym_csr_task2(C);

	int n = A->n; int m = A->m;
	int* vptr = A->vptr; int* vpos = A->vpos; double* vval = A->vval;
	int* vptr_c = C->vptr; int* vpos_c = C->vpos; double* vval_c = C->vval;
	int col_pos, row_pos;
	double* peela = malloc(m*sizeof(double));
	double* peelc = malloc(m*sizeof(double));
	char* trans = "N";

	int i;
	for ( i = 0; i < n; i++ ) {
		array_clear(peela, m);
//		array_clear(peelc, m);
		array_s2d(A, peela, i);
		mkl_cspblas_dcsrgemv(trans, &m, vval, vptr, vpos, peela, peelc);
		int k;
		for ( k = vptr_c[i]; k < vptr_c[i+1]; k++ ) {
			col_pos = vpos_c[k];
			vval_c[k] -= peelc[col_pos];
		}
	}

	free(peela); free(peelc);
}

void dgemm_sparse_upper(hbmat_t* A, hbmat_t* B, hbmat_t* C){

	/*
	 * Check if the input matrix has properly set
	 */

	if ( A->vptr == NULL )
		hyper_sym_csr_task2(A);

//	if ( B->vptr == NULL )
//		hyper_sym_csr_task1(B);

	if ( C->vptr == NULL )
		hyper_sym_csr_task2(C);


	int m = B->m; int n = B->n;
	int* vptr = B->vptr; int* vpos = B->vpos; double* vval = B->vval;
	int* vptr_c = C->vptr; int* vpos_c = C->vpos; double* vval_c = C->vval;
	int col_pos, row_pos;
	double* peelb = malloc(m*sizeof(double));
	double* peelc = malloc(m*sizeof(double));
	char* trans = "N";

	int i;
	for ( i = 0; i < m; i++ ) {
		array_clear(peelb, m);
		array_s2d(A, peelb, i);
		mkl_cspblas_dcsrgemv(trans, &m, vval, vptr, vpos, peelb, peelc);
		int k;
		for ( k = vptr_c[i]; k < vptr_c[i+1]; k++ ) {
			col_pos = vpos_c[k];
			vval_c[k] -= peelc[col_pos];
		}
	}

	free(peelb); free(peelc);
}

void dtrsm_sparse_upper(hbmat_t* A, hbmat_t* C){

	/*
	 * Check if the input matrix has properly set
	 */
//	if ( A->vptr == NULL )
//		hyper_sym_csr_task1(A);

	if ( C->vptr == NULL )
		hyper_sym_csr_task2(C);

	int n = A->n; int m = A->m;
	int* vptr = A->vptr;
	int* vpos = A->vpos;
	double* vval = A->vval;
	double* peelb = malloc(m*sizeof(double));
	double* peelc = malloc(m*sizeof(double));
	char* trans = "N";
	double alpha = 1;
	char* matdescra = "TLNC";

	int i;
	for ( i = 0; i < n; i++ ) {
		array_clear(peelb, m);
		array_s2d(C, peelb, i);
		mkl_dcsrsv(trans, &n, &alpha, matdescra, vval, vpos, vptr, vptr+1, peelb, peelc);
		array_d2s(C, peelc, i);
	}

	free(peelb); free(peelc);
}
