#ifndef __SYMFAC_H__
#define __SYMFAC_H__


#include "hb.h"


#pragma omp task in([1]A, [1]etree) out([1]entry, [1]block)
void symbolic_csr_task(int I, int J, hbmat_t *A, int b, int *etree, int *entry, hbmat_t *block);


void ereach_csr(hbmat_t *A, int r, int *etree, vector_t* sub_row);
hbmat_t* hb2hbh_sym_etree_csr(hbmat_t *A, int b, int* etree);
void ereach_csr_p(hbmat_t *A, int r, int *etree, vector_t *sub_row, vector_t *sub_val);
hbmat_t* hb2hbh_sym_etree_csr_p(hbmat_t *A, int b, int *etree);
void hyper_sym_csr_task0(int I, int J, hbmat_t *A, int b, int *etree, int *entry);
void hyper_sym_csr_task1(hbmat_t *block);
hbmat_t* hb2hbh_hyper_sym_csr(hbmat_t *A, int b, int *etree);
void hyper_sym_csr_task2(hbmat_t *block);
hbmat_t* hb2hbh_sym_etree(hbmat_t *A, int b, int* etree);
hbmat_t *hbh2hb_sym (hbmat_t *A);


//hbmat_t* ll2b(hbmatm_t *A);
//hbmatm_t* b2ll(hbmat_t *A);
//hbmatm_t* hbbm2csrbm(hbmatm_t *A);
//void llsetdiag( hbmatm_t *A );
//void m_sync(hbmatm_t *Ahb, hbmatm_t *Acsr);


#endif // __SYMFAC_H__
