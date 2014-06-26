#ifndef __CHOL_PKERNELS_H__
#define __CHOL_PKERNELS_H__


#pragma omp task inout(A[0;m*w]) priority(2)
void fac_task( int skip, int m, int w, int t, double *A); 

#pragma omp task in(A[0;m*w]) inout(B[0;m*w]) priority(1)
void upd_task(int skipa, int skipb, int m, int w, int t, double *A, double *B);


#endif // __CHOL_PKERNELS_H__
