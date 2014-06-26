#include "chol_pkernels.h"
#include "chol_utils.h"
#include "stdio.h"


void fac_task( int skip, int m, int w, int t, double *A) {
	//printf("\nfac %i %i %i %i %x\n",skip, m, w, t, A);
	double dmone = -1.00;
	double done = 1.00;
	int tm = t * m;

	A += skip;
	int r = m - skip;

	//double *Aorig = A;	

	int k;
	for ( k =0; k < w; k+=t ) {
		int b = min(w-k, t);
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
		dtrsm_("Right", "Lower", "Transpose", "Non-unit", &r, &b, &done, A, &m, B, &m);
		//print_matrix_blocked( stdout, "A", m, m, m, m, m, t, Aorig );

		double *Bj = B;
		double *Cj = B + tm; 
		int jr = r;
		int j;
		for ( j = k+t; j < w; j+=t ) {
			int jb = min(w-j, t);
			jr -= jb;
			jr = jr > 0? jr : 0;

			double *Dj = Cj + jb;
			double *Ej = Bj + jb;

			//printf("dsyrk %i %i\n",jb,b);
			dsyrk_("Lower", "No Transpose", &jb, &b, &dmone, Bj, &m, &done, Cj, &m);
			Cj += tm + t; 
	
			//printf("dgemm %i %i %i\n", jr, jb, jb);
			dgemm_("No transpose", "Transpose", &jr, &jb, &t, &dmone, Ej, &m, Bj, &m, &done, Dj, &m);
			//print_matrix_blocked( stdout, "A", m, m, m, m, m, t, Aorig );

			Bj += t;
		}

		A = B + tm;
	}
}
			      	

void upd_task(int skipa, int skipb, int m, int w, int t, double *A, double *B) {
	//printf("\nupd %i %i %i %i %i %x %x\n", skipa, skipb, m, w, t, A, B);
	double dmone = -1.00;
	double done = 1.00;
	int tm = t*m;
	
	double *Borig = B;

	A += skipb;
	B += skipb;


	int cm = m - skipb;
	int k;
	for ( k=0; k<w; k+=t ) {
		int b = min(w-k,t);

		dsyrk_("Lower", "No Transpose", &b, &w, &dmone, A, &m, &done, B, &m);
		//printf("dsyrk\n");
		//print_matrix_blocked( stdout, "A", m, m, w, w, m, t, Borig );
	
		B += b;
		double *D = A + b;
		cm -= b;

		dgemm_("No transpose", "Transpose", &cm, &b, &w, &dmone, D, &m, A, &m, &done, B, &m);
		//printf("dgemm\n");
		//print_matrix_blocked( stdout, "A", m, m, w, w, m, t, Borig );

		B += tm;
		A = D;
	}
}

