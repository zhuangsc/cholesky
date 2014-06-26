#include "chols_kernels.h"


void potrf_sparse(hbmat_t* A) {
	int n = A->n;
	int* vptr = A->vptr; int* vpos = A->vpos; 
	double* vval = A->vval;

	int J;
	for ( J=0 ; J < n; ++J ) {
		//l22 = sqrt (a22 - l12*l12_t)
		double* a22 = &(vval[vptr[J]]);
		double* l22 = a22;
		int p_a22 = vpos[vptr[J]];

		int j;
		for ( j=0; j<J; j++ ) {
			int i;
			for ( i=vptr[j]; i<vptr[j+1]; i++ ) {
				if ( vpos[i] == p_a22 ) {
					*a22 -= vval[i]*vval[i];
				}
			}
		}
		*l22 = sqrt(*a22);
		// l32 = (a32-L32*l12)/l22
		int i;
		for ( i = vptr[J]+1; i < vptr[J+1]; i++ ) {
			double* a32 = &(vval[i]);
			double* l32 = a32;
			double l31 = 0; 
			double l12 = 0;
			double l31l12 = 0;
			int p_a32 = vpos[i];

			int j;
			for ( j = 0; j < J; j++ ) {
				l31=l12=0;
				int k;
				for ( k=vptr[j]; k<vptr[j+1]; k++) {
					if (vpos[k] < p_a22)
						continue;
					if ( vpos[k] > p_a32 )
						break;
					if ( vpos[k] == p_a22 ) {
						l12 = vval[k];
					}
					if ( vpos[k] == p_a32 ) {
						l31 = vval[k]; 
						break;
					}
				}
				l31l12 += l31*l12;
			}
			*l32 = (*a32 - l31l12)/(*l22);
		}
	}
}

void dsyrk_sparse(hbmat_t* A, hbmat_t* C){
	int n = A->n; int m = A->m;
	int* vptr = A->vptr;
	int* vpos = A->vpos;
	double* vval = A->vval;
	double* peela = malloc(m*sizeof(double));
	double* peelc = malloc(m*sizeof(double));
	char* trans = "N";
	double alpha = -1; 
	double beta = 1;
	char* matdescra = "GLNC";
	hbmat_t* A_csr = hb2csr(A);

	int i;
	for ( i = 0; i < n; i++ ) {
		array_clear(peela, m);
		array_clear(peelc, m);
		array_s2d(A_csr, peela, i);
		array_s2d(C, peelc, i);
		mkl_dcscmv(trans, &m, &n, &alpha, matdescra, vval, vpos, vptr, vptr+1, peela, &beta, peelc);
		array_d2s(C, peelc, i);
	}

	free(peela); free(peelc);
	hb_free(A_csr);
}

void dgemm_sparse(hbmat_t* A, hbmat_t* B, hbmat_t* C){
	hbmat_t* B_csr = hb2csr(B);
	int m = A->m; int n = A->n;
	int* vptr = A->vptr;
	int* vpos = A->vpos;
	double* vval = A->vval;
	double* peelb = malloc(m*sizeof(double));
	double* peelc = malloc(m*sizeof(double));
	char* trans = "N";
	double alpha = -1;
	double beta = 1;
	char* matdescra = "GLNC";

	int i;
	for ( i = 0; i < n; i++ ) {
		array_clear(peelb, m);
		array_clear(peelc, m);
		array_s2d(B_csr, peelb, i);
		array_s2d(C, peelc, i);
		mkl_dcscmv(trans, &m, &n, &alpha, matdescra, vval, vpos, vptr, vptr+1, peelb, &beta, peelc);
		array_d2s(C, peelc, i);
	}

	free(peelb); free(peelc);
	hb_free(B_csr);
}

void dtrsm_sparse(hbmat_t* A, hbmat_t* C){
	hbmat_t* C_csr = hb2csr(C);
	hbmat_t* C_hb;
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
		array_clear(peelc, m);
		array_clear(peelb, m);
		array_s2d(C_csr, peelb, i);
		array_s2d(C_csr, peelc, i);
		mkl_dcscsv(trans, &n, &alpha, matdescra, vval, vpos, vptr, vptr+1, peelb, peelc);
		array_d2s(C_csr, peelc, i);
	}

	int job[6] = {0,0,0,1,1,1};
	int info;
	mkl_dcsrcsc(job,&n,C_csr->vval,C_csr->vpos,C_csr->vptr,C->vval,C->vpos,C->vptr,&info);

	free(peelb); free(peelc);
	hb_free(C_csr);
}

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
