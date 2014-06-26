#include "chols_warm.h"
#define DIM 256
#define FILL (DIM * 0.7)
void mkl_warmup(){
	srand48(time(0));
	hbmat_t *t = malloc(sizeof(hbmat_t));
	t->m = DIM; t->n = DIM;
	t->vdiag = NULL;
	int m = t->m;
	int alpha = 1; int beta = 1;
	int *vptr = t->vptr = malloc((DIM+1) * sizeof(int));
	int *vpos = t->vpos = malloc((DIM * DIM) *sizeof(int));
	double *vval = t->vval = malloc((DIM*DIM)*sizeof(double));
	vptr[0] = 0;
	int vpos_p = 0;
	puts("warm-up");
	for ( int i = 1; i <= DIM; ++i ) {
		vptr[i] = vptr[i-1] +  FILL;
		int vp = 0;
		for ( int j = vptr[i-1]; j < vptr[i]; ++j ) {
			vpos[vpos_p] = vp;
			vval[vpos_p] = drand48();
			vp++; vpos_p++;
		}
	}

	double *x = malloc(DIM*sizeof(double));
	for(int i = 0; i < DIM; ++i)
		x[i] = drand48();
	double *y = malloc(DIM*sizeof(double));
	mkl_dcsrmv("N", &m, &m, &alpha, "GLNC", vval, vpos, vptr, vptr+1, x, &beta, y);
	mkl_dcsrsv("N", &m, &alpha, "TLNC", vval, vpos, vptr, vptr+1, x, y);
	mkl_cspblas_dcsrgemv("N", &m, vval, vptr, vpos, x, y);
	
	free(x); free(y);
	free(vptr); free(vpos); free(vval);
	free(t);
}

#if 0
void print_matrix (const hbmat_t* matrix_info, int h, char* name) {
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
#endif
