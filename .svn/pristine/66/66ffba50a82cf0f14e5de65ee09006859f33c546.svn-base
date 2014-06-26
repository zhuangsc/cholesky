hbmat_t* hbh_hyper_transpose(hbmat_t *A){
	hbmat_t *B = (hbmat_t*) malloc(sizeof(hbmat_t));
	int m = A->m;
	int n = A->n;
	int elemc = A->elemc;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	hbmat_t **vval = A->vval;
	int acc = 0;

	B->m = n; B->n = m; B->elemc = elemc;
	hbmat_t **b_vval = (hbmat_t**) malloc(elemc*sizeof(hbmat_t*));
	vector_t *trans_vptr, *trans_vpos;
	trans_vptr = vector_create_size(n);
	trans_vpos = vector_create_size(elemc);
	vector_clear(trans_vptr); 
	vector_clear(trans_vpos);
	vel_t vptr_vel, vpos_vel;

	for(int j = 0; j < m; j++){
		vptr_vel.i = trans_vpos->elemc + 1;
		vector_insert(trans_vptr, vptr_vel);
		for(int i = 0; i < n; i++){
			for(int k = vptr[i]; k < vptr[i+1]; k++){
				if (j == vpos[k-1] - 1){
				 	vpos_vel.i = i + 1;
					vector_insert(trans_vpos, vpos_vel);
					b_vval[acc++] = vval[k-1];
					break;
				}
			}
		}
	}
	vptr_vel.i = trans_vpos->elemc + 1;
	vector_insert(trans_vptr, vptr_vel);
	B->vptr = vector2int(trans_vptr);
	B->vpos = vector2int(trans_vpos);
	B->vval = b_vval;

	//printf("hyper-transpose\n");
	for(int i = 0; i < B->elemc; i++){
		//((hbmat_t**)B->vval)[i] = hb_transpose(((hbmat_t**)B->vval)[i]);
		//printf("%p\t%p\n", ((hbmat_t**)B->vval)[i], b_vval[i]);
		((hbmat_t**)B->vval)[i] = hb_transpose(b_vval[i]);
		//printf("%p\t%p\n", ((hbmat_t**)B->vval)[i], b_vval[i]);
		//printf("---------------------------------\n");
	}

	return B;
}

hbmat_t* hb2hbh_sym_etree_u(hbmat_t *A, int b, int* etree){

	hbmat_t *Ab = (hbmat_t*) malloc(sizeof(hbmat_t)); 

	int m = A->m;
	int n = A->n;
	int *vptr = A->vptr;
	int *vpos = A->vpos;
	double *vval = A->vval;

	int M = ( m + b - 1 ) / b;
	int N = ( n + b - 1 ) / b;
	int num = ((1 + M) * N) / 2;  //total number of blocks in the lower triangular matrix

	Ab->m = M; Ab->n = N; Ab->elemc = 0; Ab->vdiag = NULL;
	Ab->vptr = malloc((N+1)*sizeof(int));
	Ab->vpos = malloc(num*sizeof(int));
	Ab->vval = malloc(num*sizeof(hbmat_t*));

	hbmat_t* acchb = malloc(num*sizeof(hbmat_t));
	int acc = 0 ;
	int vpos_count = 0 ;

	if ( M==0 || N==0 ) {
		fprintf( stderr, "block size %i too large\n", b);
	}

	for ( int J = 0; J < N; ++J ) { 	
		
		int jstart = J * b;
		int jc = n - jstart;
		jc = b;
		Ab->vptr[J] = acc;

		for ( int I = 0; I < J+1; ++I ) { 	

			int base_col = J * b ;
			int base_row = I * b ;
			vel_t vptr_current;
			vector_t* sub_tree = vector_create();
			vector_t* sub_vptr = vector_create();
			vector_t* sub_vpos = vector_create();
			vector_t* sub_vval = vector_create();
			vector_clear(sub_vptr); vector_clear(sub_vpos); vector_clear(sub_vval);

			for ( int j = 0 ; j < jc; j++ ) { 	// Innermost loop: column by column

				int pos_col = base_col + j ;	// Absolute column position (0 based)
				int max_row = pos_col ;		// Maximum row position in this column (0 based)
				int bborder = base_row + b - 1;
				int current_row ;
				max_row = max_row <= bborder ? max_row : bborder ;

				int min_nz = vptr[pos_col] ;
				int	max_nz = vptr[pos_col+1] ;

				vptr_current.i = sub_vpos->elemc;
				vector_insert(sub_vptr, vptr_current);
				vector_clear(sub_tree);			

				//Padding
				if (pos_col >= n){
					if(base_row+b > n){
						vptr_current.i = j;
						vector_insert(sub_vpos, vptr_current);
						vptr_current.d = 0;
						vector_insert(sub_vval, vptr_current);
					}
					continue;
				}

				for ( int current = min_nz; current < max_nz; ++current ) {
					int status ;
					current_row = vpos[current];
					vel_t vel_current = {i : current_row} ;

					if (current_row >= base_row && current_row < max_row){
						status = vector_insert_t(sub_tree, vel_current) ;
						if (!status)
							continue;
					}
					else if (current_row == max_row){
						status = vector_insert_t(sub_tree, vel_current) ; 
						break ; //current row is the boarder, break the loop
					}
					else if (current_row > max_row){
						break ;	//current row is out of the boarder, break the loop
					}

					//Traverse the elimination tree
					for (int node = etree[current_row]; node != -1 && node <= max_row; node = etree[node]){
						if (node >= base_row){
							vel_t vel_node = {i : node} ;
							vector_insert_t(sub_tree, vel_node) ;
						}
					}
				}
				
				//Sort the vector
				vector_qsorti(sub_tree);
				for(int i = 0; i < sub_tree->elemc; i++){
					vel_t vval_current = {d : 0.0};
					vel_t vpos_c;
					int fill_in = 1;

					vpos_c.i = sub_tree->elem[i].i - base_row;
					vector_insert(sub_vpos, vpos_c);
					for (int current = min_nz; current < max_nz; current++){
						current_row = vpos[current];
						if (current_row > max_row) break;
						if (current_row == sub_tree->elem[i].i){
							vval_current.d = vval[current];
							vector_insert(sub_vval, vval_current);
							fill_in = 0;
							break;
						}
					}
					if (fill_in)
						vector_insert(sub_vval, vval_current);
				}
			}
			
			vector_free(sub_tree);
			if (sub_vpos->elemc != 0){
				acchb[acc].m = jc;
				acchb[acc].n = jc;
				acchb[acc].elemc = sub_vpos->elemc;
				acchb[acc].vdiag = NULL;
				vptr_current.i = sub_vpos->elemc;
				vector_insert(sub_vptr, vptr_current);
//				vector_printi(sub_vptr);
//				vector_printi(sub_vpos);
//				vector_printd(sub_vval);
				acchb[acc].vptr = vector2int(sub_vptr);
				acchb[acc].vpos = vector2int(sub_vpos);
				acchb[acc].vval = vector2double(sub_vval);
				Ab->vpos[vpos_count] = I;
				((hbmat_t**)Ab->vval)[vpos_count] = &(acchb[acc]);
				//((hbmat_t**)Ab->vval)[vpos_count] = acchb+acc;
				vpos_count++;
				Ab->elemc++;
				acc++;
			}
			else{
				vector_free(sub_vptr);
				vector_free(sub_vpos);
				vector_free(sub_vval);
			}
		}
	}
	Ab->vptr[N] = Ab->elemc ;
	
	//TODO Check the correctness
	Ab->vpos = (int*) realloc(Ab->vpos, acc * sizeof(int));
	Ab->vval = (hbmat_t**) realloc(Ab->vval, acc * sizeof(hbmat_t*));

	return Ab;
}




hbmat_t *hbh2hb_sym (hbmat_t *A){
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


