#ifndef __FPT_H__
#define __FPT_H__


#define I_ONE	1
#define I_MONE	-1
#define F_ONE	1.0
#define F_MONE	-1.0

#ifdef SINGLE_PRECISION

typedef float fp_t;

#else

typedef double fp_t;

#endif 

extern double d_one; 
extern float s_one; 
extern double d_mone;
extern float s_mone;
extern fp_t f_one;
extern fp_t f_mone;


#endif // __FPT_H__
