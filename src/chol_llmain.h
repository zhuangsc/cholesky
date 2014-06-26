#ifndef __CHOL_LLMAIN_H__
#define __CHOL_LLMAIN_H__


#ifdef DOUBLE_PRECISION
#define CHOL_LL			dchol_ll
#else
#define CHOL_LL			schol_ll
#endif


int schol_ll(int mt, int b, int t, float **Ah);  
int dchol_ll(int mt, int b, int t, double **Ah);  


#endif // __CHOL_LLMAIN_H__
