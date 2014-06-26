#ifndef __HBCONVRT_H__
#define __HBCONVRT_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include "hb.h"
#include "vector.h"


hbmat_t* csc2csr(hbmat_t *A);
hbmat_t* hbh2hb(hbmat_t *A);
hbmat_t* hb_transpose(hbmat_t *A);
hbmat_t* hbh_hyper_transpose(hbmat_t *A);
hbmat_t *hb2hbb(hbmat_t *A, int b);
hbmat_t* hbb2csrb(hbmat_t *A);

//#pragma omp task in([1]A) out([1]entry, [1]block)
void hb2hbh_csr_task(int I, int J, hbmat_t *A, int b, int *entry, hbmat_t *block);

//#pragma omp task in([1]A) out([1]entry, [1]block)
void hb2hbh_csc_task(int I, int J, hbmat_t *A, int b, int *entry, hbmat_t *block);

hbmat_t* hb2hbh(hbmat_t *A, int b, int is_csr);



#define FILLINS 1
pthread_mutex_t mutexhb;
int *vptr_pool, *vpos_pool;
double *vval_pool;
int vptr_unit, vpos_unit, vval_unit;
int vptr_pp, vpos_pp, vval_pp;


#endif // __HBCONVRT_H__
