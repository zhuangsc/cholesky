#ifndef __MATFPRINT_H__
#define __MATFPRINT_H__


#include <stdio.h>

#include "fptype.h"
#include "hb.h"


#ifdef DOUBLE_PRECISION

#define FPRINT_DENSE2MM 		fprint_ddense2mm
#define PRINT_DENSE2MM 			print_ddense2mm
#define FPRINT_SUFF_DENSE2MM 	fprint_suff_ddense2mm

#else

#define FPRINT_DENSE2MM 		fprint_sdense2mm
#define PRINT_DENSE2MM 			print_sdense2mm
#define FPRINT_SUFF_DENSE2MM 	fprint_suff_sdense2mm

#endif


void fprint_csr2mm(const char *fname, int m, const int *vptr, const int *vpos, const double *vval, int offs);
void fprint_csc2mm(const char *fname, int m, const int *vptr, const int *vpos, const double *vval, int offs);

void fprint_sdense2mm(const char *fname, const char *name, int m, int n, const float *A);
void fprint_ddense2mm(const char *fname, const char *name, int m, int n, const double *A);

void print_sdense2mm(FILE *f, const char *name, int m, int n, const float *A);
void print_ddense2mm(FILE *f, const char *name, int m, int n, const double *A);

void fprint_suff_sdense2mm(const char *fname, int suff, const char *name, int m, int n, const float *A);  
void fprint_suff_ddense2mm(const char *fname, int suff, const char *name, int m, int n, const double *A);

void fprint_csc(const char *fname, hbmat_t *A, int offs);
void print_csc(FILE *f, hbmat_t *A, int offs);


#endif // __MATFPRINT_H__
