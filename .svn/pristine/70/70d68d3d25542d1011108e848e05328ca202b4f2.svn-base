#ifndef __BLAS_H__
#define __BLAS_H__


float sdot_(const int *, const float *, const int *, const float *, const int *);
float snrm2_(const int *, const float *, const int *);
void sgemm_(
    const char *transa, const char *transb,
    const int *m, const int *n, const int *k,
    const float *alpha,
    const float *a, const int *lda,
    const float *b, const int *ldb,
    const float *beta,
    float *c__, const int *ldc);
void scopy_(const int *n, const float *sx, const int *incx, float *sy, const int *incy);
void sscal_(const int *n, const float *sa, float *sx, const int *incx);
void saxpy_(const int *n, const float *sa, const float *sx, const int *incx, float *sy, const int *incy);
void spotrf_(const char *trian, const int *n, float *sx, const int *incx, int *info);
void ssyrk_(const char *uplo, const char *trans, const int *n, const int *k, const float *ALPHA, const float *A, const int *lda, const float *beta, float *C, const int *ldc);
void strsm_(const char *side, const char *uplo, const char *transa, const char *diag, const int *m, const int *n, const float *alpha, const float *A, const int *lda, float *B, const int *ldb);


void dgemm_(
    const char *transa, const char *transb,
    const int *m, const int *n, const int *k,
    const double *alpha,
    const double *a, const int *lda,
    const double *b, const int *ldb,
    const double *beta,
    double *c__, const int *ldc);
double ddot_(const int *, const double *, const int *, const double *, const int *);
double dnrm2_(const int *, const double *, const int *);
void dcopy_(const int *n, const double *sx, const int *incx, double *sy, const int *incy);
void dscal_(const int *n, const double *sa, double *sx, const int *incx);
void daxpy_(const int *n, const double *sa, const double *sx, const int *incx, double *sy, const int *incy);
void dpotrf_(const char *trian, const int *n, double *sx, const int *incx, int *info);
void dsyrk_(const char *uplo, const char *trans, const int *n, const int *k, const double *ALPHA, const double *A, const int *lda, const double *beta, double *C, const int *ldc);
void dtrsm_(const char *side, const char *uplo, const char *transa, const char *diag, const int *m, const int *n, const float *alpha, const float *A, const int *lda, float *B, const int *ldb);


#endif // __BLAS_H__
