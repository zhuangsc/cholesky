#ifndef __HBEXT_H__
#define __HBEXT_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h> 

#include "hb.h"
//#include "hbm.h"
#include "vector.h"


void hb_sync(hbmat_t *A, hbmat_t *B);

int hb_setdiag(hbmat_t *A);
void csr_setdiag(hbmat_t *A);

#if 0
void hbm_setdiag(hbmatm_t *A);
void csrm_setdiag(hbmatm_t *A);
#endif
void hb_markudiag(hbmat_t *A);

void one2zero(hbmat_t* in_matrix);
int* etree(hbmat_t* in_matrix);


#endif // __HBEXT_H__
