#ifndef __CHOL_NRLMAIN_H__
#define __CHOL_NRLMAIN_H__


#ifdef DOUBLE_PRECISION
#define CHOL_NRL			dchol_nrl
#else
#define CHOL_NRL			schol_nrl
#endif


int schol_nrl(int mt, int sb, int b, int t, float **Ah);
int dchol_nrl(int mt, int sb, int b, int t, double **Ah);


#endif // __CHOL_NRLMAIN_H__
