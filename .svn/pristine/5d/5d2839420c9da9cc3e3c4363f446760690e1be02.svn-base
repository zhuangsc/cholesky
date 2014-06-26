#include "hbdebug.h"

void csrb_dense_printf(FILE *f, const char *name, hbmat_t *Acsrb) {
	int m = Acsrb->m;
	int n = Acsrb->n;
	int blockc = Acsrb->elemc;
	int b = Acsrb->b;
	if ( b == 0 ) {
		fprintf( stderr, "warning: csrb_dense_printf with block size 0\n");
	}
	int *vpos = Acsrb->vpos;
	int *vptr = Acsrb->vptr;
	double **vval = Acsrb->vval;
	int *vdiag = Acsrb->vdiag;

	fprintf( f, "hbmat %s %ix%i elemc %i fill %.4f %p :\n", name, m, n, blockc, (double) 100 * blockc / ( m * n ), Acsrb );
	int i;
	for ( i = 0; i < m; i++ ) {
		printf("%i ", vptr[i]);
	}
	printf("*\n");
	if ( vdiag != NULL ) {
		for ( i = 0; i < n; i++ ) {
			printf("%i ", vdiag[i]);
		}
		printf("+\n");
	}

	for ( i = 0; i < m ; i++ ) {
		int start = vptr[i] - 1;
		int end = vptr[i+1] - 1;
		int colc = end - start;

		int c = 0;
		int j;
		for ( j = 0; j < n ; j++ ) {
			int col = vpos[start+c];
			if ( col == j + 1 && c < colc ) {
				fprintf( f, "%p ", vval[start+c]);
				++c;
			} else fprintf(f, "%p ", NULL);
		}
		fprintf( f, "\n");
	}
}


#if 0
void csrbm_dense_printf(FILE *f, const char *name, hbmatm_t *Am, int detail) {
	int m = Am->m;
	int n = Am->n;
	int blockc = Am->elemc;
	int b = Am->b;
	if ( b == 0 ) {
		fprintf( stderr, "warning: csrbm_dense_printf with block size 0\n");
	}
	//dll_t *vpos = Am->vpos;
	int *vptr = Am->vptr;
	dll3_t *vval = Am->vval;
	int *vdiag = Am->vdiag;

	fprintf( f, "hbmatm %s %ix%i elemc %i fill %.4f %p :\n", name, m, n, blockc, (double) 100 * blockc / ( m * n ), Am );

	dll3_t *vvalel = vval->next;
	int i;
	for ( i = 0; i < m ; i++ ) {
		int start = vptr[i] - 1;
		int end = vptr[i+1] - 1;
		int colc = end - start;

		int c = 0;
		int j;
		for ( j = 0; j < n ; j++ ) {
			int col = vvalel->e.i;

			if ( col == j + 1 && c < colc ) {
				if ( detail ) {
					//fprintf( f, "%p ", vval[start+c]);
				} else {
					fprintf( f, "X ");
				}
				++c;
			
				vvalel = vvalel->next;
			} else fprintf(f, "%i ", 0);

			
		}
		fprintf( f, "\n");
	}
}
#endif



void hbb_dense_printf(FILE *f, const char *name, hbmat_t *Ahbb, int struc, int detail) {
	int *vptr = Ahbb->vptr;
	int *vpos = Ahbb->vpos;
	double **vval = Ahbb->vval;
	int m = Ahbb->m;
	int n = Ahbb->n;
	int blockc = Ahbb->elemc;
	int b = Ahbb->b;
	int *vdiag = Ahbb->vdiag;

	if ( b == 0 ) { 
		fprintf( stderr, "warning: hbb_dense_printf with block size 0\n");
	}

	fprintf(f, "hbb dense %s %ix%i blockc %i fill %.4f\n", name, m, n, blockc, (double) 100 * blockc / ( m * n ) );

	int j;
	for ( j = 0 ; j < n; j++ ) {
		fprintf( f, "%i ", vptr[j]);
	}
	fprintf( f, "*\n");

	if ( vdiag ) {
		for ( j = 0 ; j < n; j++ ) {
			fprintf( f, "%i ", vdiag[j]);
		}
		fprintf( f, "+\n");
	}

	if ( struc ) {
		vector_t *cfront = vector_create();
		vector_t *celemc = vector_create();
		for ( j = 0; j < n; j++ ) {
			vel_t vel;
			vel.i = vptr[j] - 1;
			vector_insert( cfront, vel );
	
			vel.i = 0;
			vector_insert( celemc, vel );
		}

		int c = 0;
		int i;
		for ( i = 0; i < m; i++ ) {
			int j;
			for ( j = 0; j < n ; j++ ) {
				vel_t vel = vector_get( celemc, j );
				int currc = vel.i;
				int colc = vptr[j+1] - vptr[j];

				if ( currc < colc ) {
					vel_t vel = vector_get( cfront, j );
					int start = vel.i;
					int row = vpos[start];

					if ( i == row - 1 ) {
						if ( detail ) {
							fprintf(f, "%p ", vval[start]);
						} else {
							fprintf(f, "X ");
						}

						vel_t vel;
						vel.i = start+1;
						vector_insertat ( cfront, vel, j );

						vel.i = currc+1;
						vector_insertat ( celemc, vel, j );
					} else { 
						if ( detail ) {
							fprintf( f, "%p ", NULL);
						} else { 
							fprintf(f, "0 ");
						}
					}
				} else {
					if ( detail ) {
						fprintf( f, "%p ", NULL);
					} else {
						fprintf( f, "0 ");
					}
				}
			}
			printf("\n");
		} 

		vector_free(cfront);
		vector_free(celemc);
	}
}


void hb_dense_print(const char *name, hbmat_t *Ahb) {
	int *vptr = Ahb->vptr;
	int *vpos = Ahb->vpos;
	double *vval = Ahb->vval;
	int m = Ahb->m;
	int n = Ahb->n;
	int elemc = Ahb->elemc;

	printf("dense %s %ix%i \n", name, m, n);

	vector_t *cfront = vector_create();
	vector_t *celemc = vector_create();
	int j;
	for ( j = 0; j < n; j++ ) {
		vel_t vel;
		vel.i = vptr[j] - 1;
		vector_insert( cfront, vel );

		vel.i = 0;
		vector_insert( celemc, vel );
	}
		

	int c = 0;
	int i;
	for ( i = 0; i < m; i++ ) {
		int j;
		for ( j = 0; j < n ; j++ ) {
			int currc = vector_get( celemc, j ).i;
			int colc = vptr[j+1] - vptr[j];

			if ( currc < colc ) {
				int start = vector_get( cfront, j ).i;
				int row = vpos[start];

				if ( i == row - 1 ) {
					//double val = vval[start];

					printf("X ");

					vel_t vel;
					vel.i = start+1;
					vector_insertat ( cfront, vel, j );
					vel.i = currc+1;
					vector_insertat ( celemc, vel, j );
				} else printf("0 ");
			} else printf("0 ");
		}
		printf("\n");
	} 

	vector_free(cfront);
	vector_free(celemc);
}

int hb_struc_diff(hbmat_t *A, hbmat_t *B) {
	int m = A->m;
	int n = A->n;	
	int elemc = A->elemc;
	
	if ( m!=B->m || n!=B->n || elemc!=B->elemc ) {
		return 1;
	}

	int *vposA = A->vpos;;
	int *vposB = B->vpos;

	int i;
	for ( i = 0; i < m; i++ ) {
		if ( vposA[i] != vposB[i] ) {
			printf("diff at %i!\n", i);
			return 1;
		}

	}

	return 0;
}


void hb_vecat(hbmat_t *A, int i) {
	printf("%p vector at %i:\n", A, i);
	int m = A->m;
	int n = A->n;
	int elemc = A->elemc;
	int b = A->b;

	int *vptr = A->vptr;
	int *vpos = A->vpos;
	double **bval = A->vval;
	double *val = A->vval;

	int start = vptr[i-1] - 1;
	int end = vptr[i] - 1;
	int cnt = end - start;

	int c;
	for ( c=0; c<cnt; c++ ) {
		printf("%i ", vpos[start + c]);
	}
	if ( c!=0 ) {
		printf("\n");
	}

	if ( val != NULL ) {
		if ( b == 0 ) {
			for ( c=0; c<cnt; c++ ) {
				printf("%.4f ", val[start + c]);
			}
		} else { 
			for ( c=0; c<cnt; c++ ) {
				printf("%p ", bval[start + c]);
			}
		}

		if ( c!=0 ) {
			printf("\n");
		}
	}
}


void hb_elemat(hbmat_t *A, int i, int j) {
	int m = A->m;
	int n = A->n;
	int elemc = A->elemc;

	int *vptr = A->vptr;
	int *vpos = A->vpos;
	//double *vval = A->vval;

	int cstart = vptr[j-1] - 1;
	int cend = vptr[j] - 1;
	int cnt = cend - cstart;

	int c;
	for ( c = 0 ; c<cnt; c++ ) {
		int row = vpos[cstart + c];

		if ( row == i ) {
			printf("(%i,%i) found in %p\n", i, j, A);
			return ;
		}
	}

	printf("(%i,%i) not found in %p\n", i, j, A);
}

#if 0
int hb_hbm_eq(hbmat_t *A, hbmatm_t *B) {
	int elemc = A->elemc;
	if ( elemc != B->elemc ) {
		return 1;
	}	

	int m = A->m;
	int n = A->n;
	if ( m!=B->m || n!=B->n ) {
		return 1;
	}

	int *vptr = A->vptr;
	int *mvptr = B->vptr;
	int i;
	for ( i=0; i<m+1; i++ ) {
		if ( vptr[i] != mvptr[i] ) {
			return 1;
		}
	}

	dll3_t *vvall = B->vval->next;
	int *vpos = A->vpos;
	for ( i=0; i<elemc; i++ ) {
		if ( vpos[i] != vvall->e.i ) {
			return 1;
		}
		vvall = vvall->next;
	}

#if 0
	ll_t *vvall = B->vval->next;
	double **vval = A->vval;
	for ( i=0; i<elemc; i++ ) {
		if ( vval[i] != vvall->e.p ) {
			return 1;
		}
		vvall = vvall->next;
	}
#endif
	
	return 0;
}

int hb_hbm_subset(hbmat_t *Acsr, hbmatm_t *Bcsr) {
	int elemc = Acsr->elemc;
	if ( elemc > Bcsr->elemc ) {
		return 0;
	}	

	int m = Acsr->m;
	int n = Acsr->n;
	if ( m!=Bcsr->m || n!=Bcsr->n ) {
		return 0;
	}

	int *vptr = Acsr->vptr;
	int *vpos = Acsr->vpos;
	int *mvptr = Bcsr->vptr;
	dll3_t *mvval = Bcsr->vval;

	int i;
	for ( i=0; i<m; i++ ) {
		int mstarti = mvptr[i] - 1;
		int mendi = mvptr[i+1] - 1;
		dll3_t *mcol = dll3_get( mvval, mstarti );
		dll3_t *mend = dll3_get( mvval, mendi );

		int startc = vptr[i] - 1;
		int endc = vptr[i+1] - 1;
		int c;
		for (  c = startc ; c < endc ; ++c ) {
			int col = vpos[ c ];

			while ( mcol->e.i < col && mcol != mend ) {
				mcol = mcol->next;
			}

			if ( mcol->e.i != col || mcol == mend ) {
				return 0;
			}
		}
	}

	return 1;
}


void m_print(FILE *f, const char *name, hbmatm_t *A) {
	int *vptr = A->vptr;
	dll3_t *vval = A->vval;
	int *vdiag = A->vdiag;
	dll3_t **vdiagl = A->vdiagl;
	dll3_t **vptrv = A->vptrv;
	int M = A->m;
	int N = A->n;
	int elemc = A->elemc;

	fprintf(f,"mmat %s @ %p M %i N %i elemc %i fill %.4f:\n", name, A, M, N, A->elemc, (double) 100 * elemc / (M * N));

	int c;
	for ( c = 0; c <= M; c++ ) {
		fprintf(f, "%i ", vptr[c] );
	}
	fprintf(f, "*\n");

	if ( vptrv != NULL ) {
		printf("  ");
		int c;
		for ( c = 0; c <= M; c++ ) {
			fprintf(f, "%p ", vptrv[c]->e.A );
		}
		fprintf(f, "\n");
	}


	if ( vdiag!= NULL ) {
		int c;
		for ( c = 0; c < M; c++ ) {
			fprintf(f, "%i ", vdiag[c] );
		}
		fprintf(f, "+\n");
	}

	if ( vdiagl != NULL ) {
		printf("  ");
		int c;
		for ( c = 0; c < M; c++ ) {
			fprintf(f, "%p ", vdiagl[c]->e.A );
		}
		fprintf(f, "\n");
	}

	dll3_print(f, vval);

	int j;
	for ( j = 0; j<M; j++ ) {
		dll3_t *lv = vptrv[j];
		dll3_t *lend = vptrv[j+1];

		while ( lv != lend ) {
			fprintf( f, "%p ", lv->e.A );
			lv = lv->next;
		}
		
		fprintf( f, "\n" );
	}
} 


int hbm_diff(hbmatm_t *A, hbmatm_t *B) {
	int *vptr = A->vptr;
	int m = A->m;
	int elemc = A->elemc;
	dll3_t *vval = A->vval;
	int *vdiag = A->vdiag;

	int *Bvptr = B->vptr;
	int Bm = B->m;
	int Belemc = B->elemc;
	dll3_t *Bvval = B->vval;
	int *Bvdiag = B->vdiag;

	if ( m != Bm ) {
		return 1;
	}

	if ( elemc != Belemc ) {
		return 2;
	}

	int c;
	for ( c = 0; c <= m; c++ ) {
		if ( vptr[c] != Bvptr[c] ) {
			//printf("%i : %i %i\n", c, vptr[c], Bvptr[c]);
			return 3;
		}
	}

#if 0
	for ( c = 0; c < m; c++ ) {
		if ( vdiag[c] != Bvdiag[c] ) {
			return 4;
		}
	}
#endif

	if ( dll3_diff( vval, Bvval ) ) {
		return 5;
	}


	return 0;
}


int hbm_hbcsr_diff(hbmatm_t *Ahb, hbmatm_t *Acsr) {
	hbmatm_t *Aconv = hbbm2csrbm( Ahb );


	int diff = hbm_diff(Aconv, Acsr);

	hbmatm_free( Aconv );

	return diff;
}
#endif


#if 0
void print_etree(int* etree_ptr, int col){
	printf("Elimination Tree (zero-based): \n");
	for(int i = 0; i < col; i++){
		printf("node : %d \t parent : %d\n", i, etree_ptr[i]);
	}
	printf("-------------------------------\n");
}

void print_matrix(const hbmat_t* matrix_info, int h, char* name){
	double* value = (double*) matrix_info->vval;
	hbmat_t** address = matrix_info->vval;
	printf("------------------------------------\n");
	printf("\t Displaying : %s\n", name);
	printf("Rows = %d,Columns = %d,Non-zeros = %d\n",matrix_info->m,matrix_info->n,matrix_info->elemc);
	printf("Virtual Address = %p\n", matrix_info);
	for (int i = 0; i <= matrix_info->n; i++)
		printf("%d ", matrix_info->vptr[i]);
	printf("\n");
	for (int j = 0; j < matrix_info->elemc; j++)
		printf("%d ", matrix_info->vpos[j]);
	printf("\n");
	for (int j = 0; j < matrix_info->elemc; j++){
		if(h == 0)
			printf("%f ", value[j]);
		else{
			printf("vval[%d]: %p\n", j, address[j]);
			print_matrix(address[j],0,"");
		}
	}
	printf("\n");
}

void print_address(const hbmat_t* A, int h, char* name){
	printf("Matrix %s in address %p\n", name, A);
	if (h){
		for(int i = 0; i < A->elemc; i++)
			printf("sub-matrix[%d]: \t address %p\n", i, ((hbmat_t**)A->vval)[i]);
	}
}


#endif
