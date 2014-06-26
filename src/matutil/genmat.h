#ifndef __GENMAT_H__
#define __GENMAT_H__


#ifdef SINGLE_PRECISION
#define GENMAT		sgenmat
#else
#define GENMAT		dgenmat
#endif


void sgenmat(int n, int tn, int nleft, int b, float *A); 
void dgenmat(int n, int tn, int nleft, int b, double *A); 


#endif // __GENMAT_H__
