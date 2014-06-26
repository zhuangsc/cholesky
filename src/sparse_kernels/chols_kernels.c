#include "chols_kernels.h"

#include "hbconvrt.h"
#include "array.h"


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
	hbmat_t* A_csr = csc2csr(A);

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
	hbmat_t* B_csr = csc2csr(B);
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
	hbmat_t* C_csr = csc2csr(C);
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

