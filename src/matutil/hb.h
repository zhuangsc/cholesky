#ifndef __HB_H__
#define __HB_H__


#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>


#define MAT_CSC		0
#define MAT_CSR 	1
#define OFFS_F 		1
#define OFFS_C 		0


typedef struct strhbmat {
	int m, n;
	int elemc;
	int *vptr;
	int *vpos;
	void *vval;
	int *vdiag;
	int *udiagc;
	int b;
	int type;
	struct strhbmat *trans;
	struct strhbmat *orig;
	int orig_row;
	int orig_col;
	int *e_tree;

	/*
	 * The following for hyper-matrix only
	 */
	int *vptr_pool;
	int *vpos_pool;
	double *vval_pool;
	pthread_mutex_t* mtx;

} hbmat_t;


static inline void __attribute__((always_inline)) hbmat_free(hbmat_t *m) {
	if ( m->vptr != NULL ) {
		free( m->vptr );
		m->vptr = NULL;
	}

	if ( m->vpos != NULL ) {
		free( m->vpos );
		m->vpos = NULL;
	}

	if ( m->vval != NULL ) {
		free( m->vval );
		m->vval = NULL;
	}

	if ( m->vdiag != NULL ) {
		free( m->vdiag );
	}

	free(m);
}

static inline void __attribute__((always_inline)) hbmat_deepfree(hbmat_t *m) {
	int elemc = m->elemc;
	void** vval = (void**) m->vval;

	if ( vval != NULL ) {
		int c;
		for ( c = 0; c < elemc ; c++ ) {
			free( vval[c] );
		}
	}

	if ( m->vptr != NULL ) {
		free( m->vptr );
		m->vptr = NULL;
	}

	if ( m->vpos != NULL ) {
		free( m->vpos );
		m->vpos = NULL;
	}

	if ( m->vval != NULL ) {
		free( m->vval );
		m->vval = NULL;
	}

	if ( m->vdiag != NULL ) {
		free( m->vdiag );
	}

	free(m);
}

#if 0
extern int * get_sdpos(hbmat_t * A);
//extern hbmat_t* hb2hbh(hbmat_t *A, int b);
extern hbmat_t* hb2hbb(hbmat_t *A, int b);
extern hbmat_t* csc2csr(hbmat_t *A);
extern hbmat_t* hb_cp(hbmat_t *A);
#endif

void hb_print(FILE *f, const char *name, hbmat_t *A, int full);
void hb_print_CSC(char *fname, hbmat_t *A);
void hb_print_CSC2(char *fname, hbmat_t *A);
void hb_print_dense( FILE* str, char * name, hbmat_t *A, int force );
void hb_print_struc(FILE* f, const char *name, hbmat_t *Ahb);
void hbb_print_dense( FILE* str, char * name, hbmat_t *A );
int hb_diff(hbmat_t *A, hbmat_t *B);


void hb_free(hbmat_t *A);
void hbh_free(hbmat_t *A);
void hbh_free2(hbmat_t *A);


#endif // __HB_H__
