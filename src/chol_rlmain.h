#ifndef __CHOL_RLMAIN_H__
#define __CHOL_RLMAIN_H__


#ifdef DOUBLE_PRECISION
#define CHOL_RL 		dchol_rl
#else
#define CHOL_RL 		schol_rl
#endif


int schol_rl(int mt, int b, int t, float **Ah);
int dchol_rl(int mt, int b, int t, double **Ah);


#endif // __CHOL_RLMAIN_H__
