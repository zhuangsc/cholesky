#ifndef __CHOL_CHECK_H__
#define __CHOL_CHECK_H__


#include "chol_utils.h"


int chol_check(int check, int m, int mr, int ts, double **Ah, double *Rcheck);
//void qrca_zerocheck(int m, int mr, int n, int nr, int ts, double *A, char *file, int line);


#define CHOL_OCTAVE 0


#if CHOL_OCTAVE


static inline __attribute__((always_inline)) void octave_matrix_def(char *file, char *name, int m, int mr, int n, int nr, int br, int bc, double *A) {
	FILE *fstr=fopen(file,"a");

  	fprintf(fstr,"%s = [];\n",name);
	print_matrix_blocked(fstr, name, m, mr, n, nr, br, bc, A);

	fclose(fstr);
}


static inline __attribute__((always_inline)) void octave_fin(char *file, char *a,char *r)
{
  FILE *fstr=fopen(file,"a");

  fprintf(fstr,"\nt = tril(%s);\n",a);
  fprintf(fstr,"%s = t + triu(t',1);\n",a);
  fprintf(fstr,"U = chol(%s, 'lower');\n",a);
  fprintf(fstr,"norm(tril(abs(%s)-abs(U)),1)/norm(U,1)\n",r);

  fclose(fstr);
}



#else


static inline __attribute__((always_inline)) void octave_matrix_def(char *file, char *name, int m, int mr, int n, int nr, int br, int bc, double *A) {}
static inline __attribute__((always_inline)) void octave_fin(char *file, char *a,char *q) {}


#endif




#endif // __CHOL_CHECK_H__
