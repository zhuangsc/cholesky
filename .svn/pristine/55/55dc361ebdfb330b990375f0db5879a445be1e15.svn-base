#include "hbconvrt.h"

hbmat_t* hb_transpose(hbmat_t *A){
	hbmat_t *B = (hbmat_t*) malloc(sizeof(hbmat_t));
	int m = A->m;
	int n = A->n;
	int elemc = A->elemc;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	double *vval = A->vval;
	int acc = 0;

	B->m = n; B->n = m; B->elemc = elemc;
	vector_t *trans_vptr, *trans_vpos, *trans_vval;
	trans_vptr = vector_create_size(n);
	trans_vpos = vector_create_size(elemc);
	trans_vval = vector_create_size(elemc);
	vector_clear(trans_vptr); 
	vector_clear(trans_vpos);
	vector_clear(trans_vval);
	vel_t vptr_vel, vpos_vel, vval_vel;

	for(int j = 0; j < m; j++){
		vptr_vel.i = trans_vpos->elemc + 1;
		vector_insert(trans_vptr, vptr_vel);
		for(int i = 0; i < n; i++){
			for(int k = vptr[i]; k < vptr[i+1]; k++){
				if (j == vpos[k-1] - 1){
					acc++;
				 	vpos_vel.i = i + 1;
					vector_insert(trans_vpos, vpos_vel);
					vval_vel.d = vval[k-1];
					vector_insert(trans_vval, vval_vel);
					break;
				}
			}
		}
	}
	vptr_vel.i = trans_vpos->elemc + 1;
	vector_insert(trans_vptr, vptr_vel);

	B->vptr = vector2int(trans_vptr);
	B->vpos = vector2int(trans_vpos);
	B->vval = vector2double(trans_vval);
	//printf("hb_transpose\naddr_a %p \t addr_b %p\n", A, B);
	//free(A);
	return B;
}

hbmat_t *hbh2hb (hbmat_t *A){
	hbmat_t *B = (hbmat_t*) malloc(sizeof(hbmat_t));
	//Assuming the input matrix A is lower triangular
	int M = A->m; int N = A->n;
	int elemc = A->elemc;
	int* vptr = A->vptr;
	int* vpos = A->vpos;
	hbmat_t** vval = A->vval;
	
	vector_t *b_vptr, *b_vpos, *b_vval;
	b_vptr = vector_create(); b_vpos = vector_create(); b_vval = vector_create();
	vector_clear(b_vptr); vector_clear(b_vpos); vector_clear(b_vval);
	vel_t b_vptr_vel, b_vpos_vel, b_vval_vel;
	hbmat_t* sub_matrix;
	int bs = vval[0]->m; //Block size can be determined by the rows of the first sub-matrix
	int col_counter = 0;

	for(int J = 0; J < N; J++){
		sub_matrix = vval[vptr[J]]; //Fetch the first sub-matrix in this column
		int tot_col = sub_matrix->n;
		for(int j = 0; j < tot_col; j++){
			col_counter++;
			b_vptr_vel.i = b_vpos->elemc;
			vector_insert(b_vptr, b_vptr_vel);
			for(int I = vptr[J]; I < vptr[J+1]; I++){
				int c_row = vpos[I];
				int row_offset = c_row*bs;
				sub_matrix = vval[I];
				for(int jj = sub_matrix->vptr[j]; jj < sub_matrix->vptr[j+1]; jj++){
					if(1 || ((double*)sub_matrix->vval)[jj] != 0 ){
						b_vpos_vel.i = sub_matrix->vpos[jj] + row_offset;
						vector_insert(b_vpos, b_vpos_vel);
						b_vval_vel.d = ((double*)sub_matrix->vval)[jj];
						vector_insert(b_vval, b_vval_vel);
					}
				}
			}
		}
	}

	b_vptr_vel.i = b_vpos->elemc;
	vector_insert(b_vptr, b_vptr_vel);
	
//	vector_printi(b_vptr);vector_printi(b_vpos);vector_printd(b_vval);

	B->m = B->n = col_counter;
	B->elemc = b_vpos->elemc;
	B->vptr = vector2int(b_vptr);
	B->vpos = vector2int(b_vpos);
	B->vval = vector2double(b_vval);
	return B;
}

hbmat_t *csc2csr(hbmat_t *A){
	int m = A->m; int n = A->n; int elemc = A->elemc;
	int* vptr = A->vptr; int* vpos = A->vpos;
	double* vval = (double*)A->vval;
	hbmat_t* B = (hbmat_t*)malloc(sizeof(hbmat_t));
	vector_t *vptr_b, **vpos_b, **vval_b;
	vptr_b = vector_create();
	vpos_b = malloc(m*sizeof(vector_t*));
	vval_b = malloc(m*sizeof(vector_t*));
	vector_clear(vptr_b);
	vel_t b_vel;

	for(int i = 0; i < m; i++){
		vpos_b[i] = vector_create();
		vector_clear(vpos_b[i]);
		vval_b[i] = vector_create();
		vector_clear(vval_b[i]);
	}

	for (int j = 0; j < n; j++){
		for(int i = vptr[j]; i < vptr[j+1]; i++){
			b_vel.i = j;
			vector_insert(vpos_b[vpos[i]], b_vel);
			b_vel.d = vval[i];
			vector_insert(vval_b[vpos[i]], b_vel);
		}
	}

	//Update vptr_b
	b_vel.i = 0;
	vector_insert(vptr_b, b_vel);
	for(int i = 0; i < m; i++){
		b_vel = vector_get(vptr_b,i);
		b_vel.i += vpos_b[i]->elemc;
		vector_insert(vptr_b, b_vel);
	}

	for(int i = 1; i < m; i++){
		vpos_b[0] = vector_append(vpos_b[0],vpos_b[i]->elem, vpos_b[i]->elemc);
		vector_free(vpos_b[i]);
		vval_b[0] = vector_append(vval_b[0],vval_b[i]->elem, vval_b[i]->elemc);
		vector_free(vval_b[i]);
	}

	B->m = m; B->n = n; B->elemc = elemc;
	B->vptr = vector2int(vptr_b);
	B->vpos = vector2int(vpos_b[0]);
	B->vval = vector2double(vval_b[0]);
	free(vpos_b); free(vval_b);

	return B;
}




void hb2hbh_csr_task(int I, int J, hbmat_t *A, int b, int *entry, hbmat_t *block){
	int m = A->m; int n = A->n;
	int* vptr = A->vptr; int* vpos = A->vpos; 
	double* vval = A->vval;
	int brow = I*b; int erow = (I+1)*b;
	int bcol = J*b; int ecol = (J+1)*b;

	vector_t* ab_vptr = vector_create();
	vector_t* ab_vpos = vector_create();
	vector_t* ab_vval = vector_create();
	vector_clear(ab_vptr); vector_clear(ab_vpos); vector_clear(ab_vval);
	vel_t pos_val;

	*entry = 0;

	for ( int L = brow; L < erow; ++L ) {
		pos_val.i = ab_vpos->elemc; 
		vector_insert(ab_vptr, pos_val);
		
		//Padding
		if ( L >= m ) {
			if ( ecol < n )
				continue;
			pos_val.i = L - bcol;
			vector_insert(ab_vpos, pos_val);
			pos_val.d = 1;
			vector_insert(ab_vval, pos_val);
			continue;
		}

		int p_elemc = ab_vpos->elemc;
		for ( int k = vptr[L]; k < vptr[L+1]; ++k ) {
			if ( vpos[k] >= bcol && vpos[k] < ecol){
				pos_val.i = vpos[k] - bcol;
				vector_insert(ab_vpos, pos_val);
				pos_val.d = vval[k];
				vector_insert(ab_vval, pos_val);
			}
		}
	}

	pos_val.i = ab_vpos->elemc;
	vector_insert(ab_vptr, pos_val);

	if ( ab_vpos->elemc ){
		block->m = b; block->n = b; block->elemc = ab_vpos->elemc;
		block->vdiag = NULL;
		block->vptr = vector2int(ab_vptr); 
		block->vpos = vector2int(ab_vpos);
		block->vval = vector2double(ab_vval);
		*entry = 1;
	}
}

void hb2hbh_csc_task(int I, int J, hbmat_t *A, int b, int *entry, hbmat_t *block){
	int m = A->m ;int n = A->n;
	int* vptr = A->vptr; int* vpos = A->vpos; 
	double* vval = A->vval;
	int brow = I*b; int erow = (I+1)*b;
	int bcol = J*b; int ecol = (J+1)*b;

	vector_t* ab_vptr = vector_create();
	vector_t* ab_vpos = vector_create();
	vector_t* ab_vval = vector_create();
	vector_clear(ab_vptr); vector_clear(ab_vpos); vector_clear(ab_vval);
	vel_t pos_val;

	*entry = 0;

	for ( int L = bcol; L < ecol; ++L ) {
		pos_val.i = ab_vpos->elemc; 
		vector_insert(ab_vptr, pos_val);
		
		//Padding
		if ( L >= n ) {
			if ( erow < m )
				continue;
			pos_val.i = L - brow;
			vector_insert(ab_vpos, pos_val);
			pos_val.d = 1;
			vector_insert(ab_vval, pos_val);
			continue;
		}

		int p_elemc = ab_vpos->elemc;
		for ( int k = vptr[L]; k < vptr[L+1]; ++k ) {
			if ( vpos[k] >= brow && vpos[k] < erow){
				pos_val.i = vpos[k] - brow;
				vector_insert(ab_vpos, pos_val);
				pos_val.d = vval[k];
				vector_insert(ab_vval, pos_val);
			}
		}
	}

	pos_val.i = ab_vpos->elemc;
	vector_insert(ab_vptr, pos_val);

	if ( ab_vpos->elemc ){
		block->m = b; block->n = b; block->elemc = ab_vpos->elemc;
		block->vdiag = NULL;
		block->vptr = vector2int(ab_vptr); 
		block->vpos = vector2int(ab_vpos);
		block->vval = vector2double(ab_vval);
		*entry = 1;
	}
}

hbmat_t* hb2hbh(hbmat_t *A, int b, int is_csr){

	int m = A->m; int n = A->n; int elemc = A->elemc;
	int *vptr = A->vptr; int *vpos = A->vpos; double* vval = A->vval;
	int M = ( m + b - 1 ) / b;
	int N = ( n + b - 1 ) / b;
	int num = M * N;

	hbmat_t* hyper = malloc(sizeof(hbmat_t));
	hyper->m = M; hyper->n = N; hyper->vdiag = NULL;
	hyper->vval = malloc(num * sizeof(hbmat_t*));
	hbmat_t** hbmat_array = malloc(num * sizeof(hbmat_t*));
	int* hentry = malloc(num * sizeof(int));

	vector_t* ab_vptr = vector_create(); 
	vector_t* ab_vpos = vector_create();
	vector_clear(ab_vptr); vector_clear(ab_vpos);
	vel_t pos_val;

	if ( M==0 || N==0 ) {
		fprintf( stderr, "block size %i too large\n", b);
	}

	for(int i = 0; i < num; ++i)
		hbmat_array[i] = malloc(sizeof(hbmat_t));

	int acc = 0;
	int I, J;
	if (is_csr){
		for ( I = 0; I < M; ++I ) {
			for ( J = 0; J < N; ++J ) {
				hb2hbh_csr_task(I, J, A, b, &(hentry[acc]), hbmat_array[acc]);
				++acc;
			}
		}
	}else{
		for ( J = 0; J < N; ++J ) {
			for ( I = 0; I < M; ++I ) {
				hb2hbh_csc_task(I, J, A, b, &(hentry[acc]), hbmat_array[acc]);
				++acc;
			}
		}
	}

//#pragma omp taskwait

	acc = 0;
	int acc0 = 0;
	if ( is_csr ){
		for ( I = 0; I < M; ++I ) {
			pos_val.i = ab_vpos->elemc;
			vector_insert(ab_vptr, pos_val);
			for ( J = 0; J < N; ++J ) {
				if ( hentry[acc] ) {
					pos_val.i = J;
					vector_insert(ab_vpos, pos_val);
					((hbmat_t**)hyper->vval)[acc0] = hbmat_array[acc];
					++acc;
					++acc0;
				} else {
					free(hbmat_array[acc]);
					++acc;
				}
			}
		}
	} else {
		for ( J = 0; J < N; ++J ) {
			pos_val.i = ab_vpos->elemc;
			vector_insert(ab_vptr, pos_val);
			for ( I = 0; I < M; ++I ) {
				if ( hentry[acc] ) {
					pos_val.i = I;
					vector_insert(ab_vpos, pos_val);
					((hbmat_t**)hyper->vval)[acc0] = hbmat_array[acc];
					++acc;
					++acc0;
				} else {
					free(hbmat_array[acc]);
					++acc;
				}
			}
		}
	}

	pos_val.i = ab_vpos->elemc;
	vector_insert(ab_vptr, pos_val);
	hyper->elemc = ab_vpos->elemc;
//	vector_printi(ab_vptr);
//	vector_printi(ab_vpos);
	hyper->vptr = vector2int(ab_vptr);
	hyper->vpos = vector2int(ab_vpos);

	set_diag(hyper);

	return hyper;
}


void set_diag(hbmat_t *A){
	int m = A->m; int n = A->n;
	int *vdiag = calloc(n, sizeof(int));
	A->vdiag = vdiag;
	int *vptr = A->vptr; int *vpos = A->vpos;
	for (int j = 0; j < n; ++j) {
		for (int k = vptr[j]; k < vptr[j+1]; ++k) {
			if (vpos[k] == j)
				vdiag[j] = k;
		}
	}
}

