#include "matfprint.h"

#include <stdio.h>
#include <string.h>
#include <math.h>


#if 0
fprint_hb2mm(const char *fname, const hbmat_t *Ahb) {
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	double *vval = A->vval;
	int n = A->n;
}
#endif


void fprint_csr2mm(const char *fname, int m, const int *vptr, const int *vpos, const double *vval, int offs) {
	int roffs = 1 - offs;
	FILE *f = fopen(fname, "w");
	
	int i;
	for ( i=0; i<m; ++i ) {
		int rowb = vptr[i] - offs;
		int rowe = vptr[i+1] - offs;
		int j;
		for ( j=rowb; j<rowe; ++j ) {
			if ( vval == NULL ) {
				fprintf(f, "%i\t%i\t1\n", i+1, vpos[j] + roffs);
			} else {
				fprintf(f, "%i\t%i\t%.16e\n", i+1, vpos[j] + roffs, vval[j]);
			}
		}
	}

	fclose(f);
}


/* can be done easier if matrix is symmetric */
void fprint_csc2mm(const char *fname, int m, const int *vptr, const int *vpos, const double *vval, int offs) {
	int roffs = 1 - offs;
	FILE *f = fopen(fname, "w");

	int *front = malloc(m * sizeof(int));
	memcpy(front, vptr, m * sizeof(int));
	
	int done = 0;
	int i = 0;
	while ( done != m ) {
		int j;
		for ( j=0; j<m; ++j ) {
			if ( front[j] >= 0 ) {
				if ( front[j] == vptr[j+1] ) {
					front[j] = -1;
					++done;
				} else {
					int x = front[j] - offs;
					int r = vpos[x] - offs;
					if ( r == i ) {
						fprintf(f, "%i\t%i\t%.16e\n", r+1, j+1, vval[x]);
						++front[j];  
					}
				}
			}
		}
		++i;
	}

	free(front);

	fclose(f);
}


void print_sdense2mm(FILE *f, const char *name, int m, int n, const float *A) {
	fprintf(f, "# name: %s\n", name);
	fprintf(f, "# type: matrix\n");
	fprintf(f, "# rows: %i\n", m);
	fprintf(f, "# columns: %i\n", n);

	int i;
	for ( i=0; i<m; ++i ) {
		int j;
		for ( j=0; j<n; ++j ) {
				fprintf(f, "%.16e ", A[j*m+i]);
		}
		fprintf(f, "\n");
	}
}


void fprint_sdense2mm(const char *fname, const char *name, int m, int n, const float *A) {
	FILE *f = fopen(fname, "w");
	if ( f == NULL ) {
		fprintf(stderr, "err: cannot open %s for writing\n", fname);
	}

	print_sdense2mm(f, name, m, n, A);

	fclose(f);
}


void print_ddense2mm(FILE *f, const char *name, int m, int n, const double *A) {
	printf("warning: writing obj %s\n", name);

	fprintf(f, "# name: %s\n", name);
	fprintf(f, "# type: matrix\n");
	fprintf(f, "# rows: %i\n", m);
	fprintf(f, "# columns: %i\n", n);

	int i;
	for ( i=0; i<m; ++i ) {
		int j;
		for ( j=0; j<n; ++j ) {
				fprintf(f, "%.16e ", A[j*m+i]);
		}
		fprintf(f, "\n");
	}
}


void fprint_ddense2mm(const char *fname, const char *name, int m, int n, const double *A) {
	FILE *f = fopen(fname, "w");
	if ( f == NULL ) {
		fprintf(stderr, "err: cannot open %s for writing\n", fname);
	}

	print_ddense2mm(f, name, m, n, A);

	fclose(f);
}


void fprint_suff_sdense2mm(const char *fname, int suff, const char *name, int m, int n, const float *A) {
		char buf1[128];
		char buf2[128];
		sprintf(buf1, "%s_%i.mm", name, suff);
		sprintf(buf2, "%s_%i", name, suff);

		fprint_sdense2mm(buf1, buf2, m, n, A);
}

void fprint_suff_ddense2mm(const char *fname, int suff, const char *name, int m, int n, const double *A) {
		char buf1[128];
		char buf2[128];
		sprintf(buf1, "%s_%i.mm", name, suff);
		sprintf(buf2, "%s_%i", name, suff);

		fprint_ddense2mm(buf1, buf2, m, n, A);
}


void print_csc(FILE *f, hbmat_t *A, int offs) {
	int roffs = 1 - offs;
	int m = A->m;
	int n = A->n;
	int elemc = A->elemc;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	double *vval = A->vval;

	int nanflag = 0;
	int j;
	for ( j=0; j<n; ++j ) {
		int b = vptr[j] - offs; 
		int e = vptr[j+1] - offs;
        int i;
        for ( i=b; i<e; ++i ) {
				double val = vval[i];
				nanflag |= isnan(val);
				
				fprintf(f, "%i %i %.16f\n", vpos[i]+roffs, j+1, vval[i]);
        }
	}

	if ( nanflag ) {
		fprintf(stderr, "warn: found NaN\n");
	}
}


void fprint_csc(const char *fname, hbmat_t *A, int offs) {
	FILE *f = fopen(fname, "w");
	if ( f == NULL ) {
		fprintf(stderr, "error: cannot open %s for writing\n", fname);
	}

	print_csc(f, A, offs);

	fclose(f);
}


