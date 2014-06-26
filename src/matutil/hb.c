#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>


#include "hb.h"
#include "vector.h"


#if 0
int ProdSparseMatrixVector (SparseMatrix spr, double *vec, double *res)
	{
		int i, j;
		double aux;

		for (i=0; i<spr.dim1; i++)
			{
				aux = 0.0;
				for (j=spr.vptr[i]; j<spr.vptr[i+1]; j++)
					aux += spr.vval[j] * vec[spr.vpos[j]];
				res[i] += aux;
			}
	}

int ProdSymSparseMatrixVector (SparseMatrix spr, double *vec, double *res)
	{
		int i, j, k;
		double aux, val;

		for (i=0; i<spr.dim1; i++)
			{
				aux = 0.0;
				for (j=spr.vptr[i]; j<spr.vptr[i+1]; j++) {
					k = spr.vpos[j]; val = spr.vval[j];
					aux += val * vec[k];
					if (k != i) res[k] += (val * vec[i]);
				}
				res[i] += aux;
			}
	}
#endif


#if 0
void CreateSparseMatrixHB (char *nameFile, ptr_SparseMatrix spr, int FtoC)
	{
		FILE *file;
		char string[length1], *s = NULL;
		int i, j, k = 0, shft = (FtoC)?-1:0;
		int *vptr = NULL, *vpos = NULL;
		double *vval = NULL; 
		int lines[5], dim[4], formats[10];

		file = hb_open (nameFile, "r");
		ReadStringFile (file, string, length1);
		ReadStringFile (file, string, length1);
		GetIntsFromString (string, lines, 5, 14, 0); 
		ReadStringFile (file, string, length1);
		GetIntsFromString ((string+14), dim, 4, 14, 0);

		CreateSparseMatrix (spr, dim[0], dim[1], dim[2], 0);
		vptr = spr->vptr; vpos = spr->vpos; vval = spr->vval; 

		ReadStringFile (file, string, length1);
		GetFormatsFromString (string, formats, 2, 16);
		GetFormatsFromString ((string+32), (formats+4), 1+(lines[4] > 0), 20);

		if (lines[4] > 0) ReadStringFile (file, string, length1);

		j = 0;
		for (i = 0; i < lines[1]; i++)
			{
				ReadStringFile (file, string, length2); 
				k = ((dim[0] + 1) - j);
				if (k > formats[0]) k = formats[0];
				GetIntsFromString (string, (vptr+j), k, formats[1], shft);
				j+=formats[0];
			}

		j = 0;
		for (i = 0; i < lines[2]; i++)
			{
				ReadStringFile (file, string, length2); 
				k = (dim[2] - j);
				if (k > formats[2]) k = formats[2];
				GetIntsFromString (string, (vpos+j), k, formats[3], shft);
				j+=formats[2];
			}

		int elemc = 0;
		j = 0;
		for (i = 0; i < lines[3]; i++)
			{
				ReadStringFile (file, string, length2); 
				k = (dim[2] - j);
				if (k > formats[4]) k = formats[4];
				GetDoublesFromString (string, (vval+j), k, formats[5]);
				elemc += k;
				
				j+=formats[4];
			}

		fclose (file);

	spr->elemc = elemc;
	printf("elemc %i\n", elemc);
}
#endif


int * get_sdpos(hbmat_t * A) {
	int n = A->n;
	int *vptr = A->vptr;
	int *vpos = A->vpos;

	int *sdpos = (int*) calloc(n, sizeof(int));

	int j;
	for ( j=0; j<n; j++ ) {
		int jj = j + 1;

		int start = vptr[j];
		int elemc = vptr[j+1] - start;
		--start;
		
		int pos = 0;
		int i;
		while ( i<elemc && pos<jj ) {
			pos = vpos[start++];
			++i;
		}

		if ( pos > jj ) {
			sdpos[j] = pos;
		}
	}

	return sdpos;
}


void hb_print_CSC(char *fname, hbmat_t *A) {
        int m = A->m;
        int n = A->n;
        int elemc = A->elemc;
        int *vptr = A->vptr;
        int *vpos = A->vpos;
        double *vval = A->vval;

        FILE *f = fopen(fname, "w");
	if ( f == NULL ) {
		fprintf(stderr, "error: cannot open %s for writing\n", fname);
	}

        
        fprintf(f, "# name: r\n");
        fprintf(f, "# type: matrix\n");
        fprintf(f, "# rows: %i\n", elemc);
        fprintf(f, "# columns: 1\n");
        int j;
        for ( j=0; j<elemc; j++ ) {
                fprintf(f, " %i\n", vpos[j]);
        }

        fprintf(f, "\n# name: c\n");
        fprintf(f, "# type: matrix\n");
        fprintf(f, "# rows: %i\n", elemc);
        fprintf(f, "# columns: 1\n");
        for ( j=0; j<n; j++ ) {
                int vc = vptr[j+1] - vptr[j]; 
                int jj;
                for ( jj=0; jj<vc; jj++) {
                        fprintf(f, " %i\n", j+1);
                }
        }

        fprintf(f, "\n# name: v\n");
        fprintf(f, "# type: matrix\n");
        fprintf(f, "# rows: %i\n", elemc);
        fprintf(f, "# columns: 1\n");
        for ( j=0; j<elemc; j++ ) {
                fprintf(f, " %.16f\n", vval[j]);
        }

	fclose(f);
}

void hb_print_CSC2(char *fname, hbmat_t *A) {
	int n = A->n;
	int m = A->m;
	int elemc = A->elemc;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	double *vval = A->vval;

	FILE *f = fopen(fname, "w");
	if ( f == NULL ) {
		fprintf(stderr, "error: cannot open %s for writing\n", fname);
	}

//	fprintf(f, "# name: %s\n", vname);
//	fprintf(f, "# type: sparse matrix\n");
//	fprintf(f, "# nnz: %d\n", elemc);
//	fprintf(f, "# rows: %d\n# columns: %d\n", m, n);
	int j;
	for ( j=0; j<n; ++j ) {
		int b = vptr[j]; 
		int e = vptr[j+1];
        int i;
        for ( i=b; i<e; ++i ) {
                //fprintf(f, "%i\t%i\t%.16f\n", vpos[i], j+1, vval[i]);
				fprintf(f, "%i %i %.16f\n", vpos[i]+1, j+1, vval[i]);
        }
	}

	fclose(f);
}

void hb_print_dense( FILE* str, char * name, hbmat_t *A, int force ) {
	fprintf(str, "%s = [];\n", name);

	int * vptr = A->vptr;
	int * vpos = A->vpos;
	double * vval = A->vval;
	int m = A->m;
	int n = A->n;

	int j;
	for ( j = 0; j < n; j++ ) {
		int start = vptr[j];
		int elemc = vptr[j+1] - start;
		--start;
	
		int c = 0;
		int i;
		for ( i = 0; i < m ; i++ ) {
			double b = 0.0;
			int pos = start + c;
			if ( vpos[pos] == i+1 && c < elemc ) {
				if ( vval == NULL || force) {
					b = 1.0;
				} else {
					b = vval[pos];
					if ( b == 0 ) { 
						fprintf(stderr, "warning: found 0 in non-zero entry!\n");
					}
				}
				++c;
			}
			fprintf(str, "%s(%d,%d) = %.14f;\n", name, i+1, j+1, b);
		}
	}
}


void hb_print_struc(FILE* f, const char *name, hbmat_t *Ahb) {
	int *vptr = Ahb->vptr;
	int *vpos = Ahb->vpos;
	double *vval = Ahb->vval;
	int m = Ahb->m;
	int n = Ahb->n;
	int elemc = Ahb->elemc;
	int b = Ahb->b;

	fprintf(f, "dense %s %ix%i b %i\n", name, m, n, b);

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

					fprintf(f, "X ");

					vel_t vel;
					vel.i = start+1;
					vector_insertat ( cfront, vel, j );
					vel.i = currc+1;
					vector_insertat ( celemc, vel, j );
				} else fprintf(f, "  ");
			} else fprintf(f, "  ");
		}
		fprintf(f, "\n");
	} 

	vector_free(cfront);
	vector_free(celemc);
}


void hbb_print_dense( FILE* str, char * name, hbmat_t *A ) {
	fprintf(str, "%s = [];\n", name);

	int * vptr = A->vptr;
	int * vpos = A->vpos;
	double **vval = A->vval;
	int M = A->m;
	int N = A->n;
	int b = A->b;

	//fprintf(stderr, "warning: hbb_print_dense: printing matrix %ix%i b %i\n", M, N, b);

	double *debug = vval[0];

	int J;
	for ( J = 0; J < N; J++ ) {
		int start = vptr[J];
		int blockc = vptr[J+1] - start;
		--start;
	
		int C = 0;
		int I;
		for ( I = 0; I < M ; I++ ) {
			int pos = start + C;
			if ( vpos[pos] == I+1 && C < blockc ) {
				double *B = vval[pos];

				int j;
				for ( j = 0; j < b; j++ )  {
					int i;
					for ( i = 0 ; i < b; i++ ) {
						int r = I * b + i + 1;
						int c = J * b + j + 1;
			
						fprintf(str, "%s(%d,%d) = %.14f;\n", name, r, c, B[j*b+i]);
					}
				}
				
				++C;
			}
		}
	}
}


void hb_print(FILE *f, const char *name, hbmat_t *A, int full) {
	fprintf(f, "hbmat %s %p\n", name, A);

	int M = A->m;
	int N = A->n;
	int b = A->b;
	fprintf(f, "M %i N %i elemc %i fill %.4f\n", M, N, A->elemc, (double) 100 * A->elemc / (M * N));
	
	int *vptr = A->vptr;
	int *vpos = A->vpos;

	int j;
	for ( j=0; j<=N; j++ ) {
		fprintf(f, "%i ", vptr[j]);
	}
	printf("*\n");

	int *vdiag = A->vdiag;
	if ( vdiag != NULL ) {
		int j;
		for ( j=0; j<N; j++ ) {
			fprintf(f, "%i ", vdiag[j]);
		}
		fprintf(f, "+\n");
	}
	
	if ( full ) {
		for ( j=0; j<N; j++ ) {
			int start = vptr[j];
			int elemc = vptr[j+1] - start;
			--start; 

			int i;
			for ( i=0; i< elemc; i++ ) {
				printf("%i ", vpos[start+i] );
			}
			printf(" )\n");
		}

		if ( b == 0 ) {
			double *vval = A->vval;
			for ( j=0; j<N; j++ ) {
				int start = vptr[j];
				int elemc = vptr[j+1] - start;
				--start; 

				int i;
				for ( i=0; i< elemc; i++ ) {
					printf("%.4f ", vval[start+i] );
				}
				printf("\n");
			}
			printf("\n");
		}
	}
}


int hb_diff(hbmat_t *A, hbmat_t *B) {
	int m = A->m;
	int n = A->n;

	if ( m != B->m ) {
		return 1;
	}

	if ( n != B->n ) {
		return 2;
	}

	int elemc = A->elemc;
	if ( elemc != B->elemc ) {
		return 3;
	}

	int *vptr = A->vptr;
	int *Bvptr = B->vptr;
	int j;
	for ( j=0; j<=n; j++ ) {
		if ( vptr[j] != Bvptr[j] ) {
			return 4;
		}
	}

	int *vval = A->vval;
	int *Bvval = B->vval;
	int c;
	for ( c=0; c<elemc; c++ ) {
		if ( vval[c] != Bvval[c] ) {
			return 5;
		}
	}

	return 0;
}


hbmat_t* hb_cp(hbmat_t *A) {
	int m = A->m;
	int n = A->n;
	int elemc = A->elemc;

	hbmat_t *cp = malloc( sizeof(hbmat_t) );
	cp->m = m;
	cp->n = n;
	cp->elemc = elemc;
	cp->b = -1;

	cp->vptr = malloc( sizeof(int) * n );
       	memcpy(cp->vptr, A->vptr, n * sizeof(int));
	cp->vval = malloc( sizeof(double) * elemc );
       	memcpy(cp->vval, A->vval, elemc * sizeof(double));
	cp->vpos = malloc( sizeof(int) * elemc );
       	memcpy(cp->vval, A->vpos, elemc * sizeof(double));

	cp->vdiag = NULL;

	return cp;
}


void hb_free(hbmat_t *A){
	free(A->vptr); free(A->vpos);
	free(A->vval);
	free(A);
}

void hbh_free(hbmat_t *A){
	int elemc = A->elemc;
	for(int i = 0; i < elemc; i++){
		free(((hbmat_t**)A->vval)[i]->vptr);
		free(((hbmat_t**)A->vval)[i]->vpos);
		free(((hbmat_t**)A->vval)[i]->vval);
	}
	free(((hbmat_t**)A->vval)[0]);
	hb_free(A);
}

void hbh_free2(hbmat_t *A){

	pthread_mutex_destroy(A->mtx);
	int elemc = A->elemc;
	free(A->e_tree);
	free(((hbmat_t**)A->vval)[0]);
	free(A->vptr_pool);
	free(A->vpos_pool);
	free(A->vval_pool);
	hb_free(A);

}


