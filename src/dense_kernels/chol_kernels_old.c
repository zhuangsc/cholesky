#include "chol_kernels.h"

#include <stdio.h>

#include "fptype.h"
#include "fpblas.h"
#include "blas.h"


void GEMM_TASK( int m, int t, fp_t *A, fp_t *B, fp_t *C) {
#if 1
	BLAS_gemm(MAT_NOTRANSP, MAT_TRANSP, m, m, m, f_mone, A, m, B, m, f_one, C, m);
#else
	int tm = t*m;
	int k;
	for ( k=0; k<m; k+=t ) {
		int b = min(m-k,t);

			dgemm_("No transpose", "Transpose", &m, &t, &m, &mone, A, &m, B, &m, &one, C, &m);

		C += tm;
		B += t;
	}
#endif
}
			      	

void SYRK_TASK(int w, fp_t *A, fp_t *C) {
	BLAS_syrk(MAT_LTRIAN, MAT_NOTRANSP, w, w, f_mone, A, w, f_one, C, w);
}


void POTRF_TASK( int m, int t, fp_t *A) {
#if 1
	int info;
	BLAS_potrf(MAT_LTRIAN, m, A, m, info);
	if ( info!=0 ) {
		fprintf(stderr,"error in POTRF (%i)\n",info);
	}
#else
	int tm = t * m;

	int r = m; 

	int k;
	for ( k =0; k < m; k+=t ) {
		int b = min(m-k, t);
		double *B = A + b;
		r -= b;

		//printf("dpotrf b %i\n", b);
		int info;
		dpotrf_("Lower", &b, A, &m, &info);
		if (info!=0) {
			fprintf(stderr,"error: DPOTRF (%i)\n", info);
		}

		r = r>0 ? r : 0;
		//printf("dtrsm %i %i\n", r, b);
		dtrsm_("Right", "Lower", "Transpose", "Non-unit", &r, &b, &d_one, A, &m, B, &m);
		//print_matrix_blocked( stdout, "A", m, m, m, m, m, t, Aorig );

		double *Bj = B;
		double *Cj = B + tm; 
		int jr = r;
		int j;
		for ( j = k+t; j < m; j+=t ) {
			int jb = min(m-j, t);
			jr -= jb;
			jr = jr > 0? jr : 0;

			double *Dj = Cj + jb;
			double *Ej = Bj + jb;

			//printf("dsyrk %i %i\n",jb,b);
			dsyrk_("Lower", "No Transpose", &jb, &b, &d_mone, Bj, &m, &d_one, Cj, &m);
			Cj += tm + t; 
	
			//printf("dgemm %i %i %i\n", jr, jb, jb);
			dgemm_("No transpose", "Transpose", &jr, &jb, &t, &d_mone, Ej, &m, Bj, &m, &d_one, Dj, &m);
			//print_matrix_blocked( stdout, "A", m, m, m, m, m, t, Aorig );

			Bj += t;
		}

		A = B + tm;
	}
#endif
}
			      	

void TRSM_TASK( int m, int t, fp_t *A, fp_t *B) {
#if 1
	BLAS_trsm(MAT_RIGHT, MAT_LTRIAN, MAT_TRANSP, MAT_DIAG_NUNIT, m, m, f_one, A, m, B, m);
#else
	int tm = t*m;
	int km = m;
	int k; 
	for ( k=0; k<m; k+=t ) {
		int b = min(m-k,t);
		
		dtrsm_("Right", "Lower", "Transpose", "Non-unit", &m, &b, &d_one, A, &m, B, &m);

		double *C = B + tm;
		km -= t;
		A += t;

		dgemm_("No Transpose", "Transpose", &m, &km, &t, &d_mone, B, &m, A, &m, &d_one, C, &m);

		B = C ;
		A += tm;
	}
#endif
}
