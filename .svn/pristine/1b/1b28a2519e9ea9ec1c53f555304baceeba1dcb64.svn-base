#ifndef __HBDEBUG_H__
#define __HBDEBUG_H__


#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>


#include "hb.h"
//#include "hbm.h"
#include "hbconvrt.h"
//#include "dlinkdlist3.h"
#include "vector.h"


//void hbm_print(FILE *f, const char *name, hbmatm_t *A); 
void hb_dense_print(const char *name, hbmat_t *Ahb);
int hb_struc_diff(hbmat_t *A, hbmat_t *B);

void csrb_dense_printf(FILE *f, const char *name, hbmat_t *Acsrb);
void hbb_dense_printf(FILE *f, const char *name, hbmat_t *Ahbb, int struc, int detail);
//void csrbm_dense_printf(FILE *f, const char *name, hbmatm_t *Am, int detail);

void hb_elemat(hbmat_t *A, int i, int j);
void hb_vecat(hbmat_t *A, int i);

//int hb_hbm_eq(hbmat_t *A, hbmatm_t *B);
//int hb_hbm_subset(hbmat_t *Acsr, hbmatm_t *Bcsr);
//void m_print(FILE *f, const char *name, hbmatm_t *A);

void print_etree(int* etree_ptr, int col);
void print_matrix(const hbmat_t* matrix_info, int h, char* name);
void print_address(const hbmat_t* A, int h, char* name);


#endif // __HBDEBUG_H__

