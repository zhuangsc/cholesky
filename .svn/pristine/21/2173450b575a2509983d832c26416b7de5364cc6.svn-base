#ifndef __FPBLAS_H__
#define __FPBLAS_H__


#include "config.h"


#if USE_MKL

#include "mkl.h"

#define ROW_MAJOR 			CblasRowMajor
#define COL_MAJOR			CblasColMajor
#define MAT_TRANSP			CblasTrans
#define MAT_NOTRANSP		CblasNoTrans
#define MAT_UTRIAN			CblasUpperTriang
#define MAT_LTRIAN			CblasLowerTriang

#ifdef SINGLE_PRECISION

#define BLAS_cp(n, dx, incx, dy, incy) 									cblas_scopy(n, dx, incx, dy, incy)
#define BLAS_gemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) 			cblas_sgemm(COL_MAJOR, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
#define BLAS_dot(n, dx, incx, dy, incy) 								cblas_sdot(n, dx, incx, dy, incy) 
#define BLAS_axpy(n, da, dx, incx, dy, incy) 								cblas_saxpy(n, da, dx, incx, dy, incy)
#define BLAS_scal(n, da, dx, incx)									cblas_sscal(n, da, dx, incx)
#define BLAS_nrm2(n, x, incx)										cblas_snrm2(n, x, incx)

#else

#define BLAS_cp(n, dx, incx, dy, incy) 									cblas_dcopy(n, dx, incx, dy, incy)
#define BLAS_gemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) 			cblas_dgemm(COL_MAJOR, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
#define BLAS_dot(n, dx, incx, dy, incy) 								cblas_ddot(n, dx, incx, dy, incy) 
#define BLAS_axpy(n, da, dx, incx, dy, incy) 								cblas_daxpy(n, da, dx, incx, dy, incy)
#define BLAS_scal(n, da, dx, incx)									cblas_dscal(n, da, dx, incx)
#define BLAS_nrm2(n, x, incx)										cblas_dnrm2(n, x, incx)

#endif 


#else /* NOT MKL */


#include "blas.h"

#define ROW_MAJOR 		
#define COL_MAJOR		
#define MAT_TRANSP			"T"
#define MAT_NOTRANSP		"N"
#define MAT_UTRIAN			"U"
#define MAT_LTRIAN			"L"
#define MAT_DIAG_UNIT		"U"
#define MAT_DIAG_NUNIT		"N"
#define MAT_LEFT			"L"
#define MAT_RIGHT			"R"

#ifdef SINGLE_PRECISION

#define BLAS_cp(n, dx, incx, dy, incy) 														scopy_(&n, dx, &incx, dy, &incy)
#define BLAS_gemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) 			sgemm_(transa, transb, &(m), &(n), &(k), &(alpha), A, &(lda), B, &(ldb), &(beta), C, &(ldc))
#define BLAS_dot(n, dx, incx, dy, incy) 													sdot_(&(n), dx, &(incx), dy, &(incy)) 
#define BLAS_axpy(n, da, dx, incx, dy, incy) 												saxpy_(&n, &da, dx, &incx, dy, &incy)
#define BLAS_scal(n, da, dx, incx)															sscal_(&n, &da, dx, &incx)
#define BLAS_nrm2(n, x, incx)																snrm2_(&n, x, &incx)
#define BLAS_syrk(trian, trans, n, k, alpha, A, lda, beta, C, ldc)							ssyrk_(trian, trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc)		
#define BLAS_potrf(trian, n, A, lda, info)													spotrf_(trian, &n, A, &lda, &info)
#define BLAS_trsm(side, trian, trans, diag, m, n, alpha,A, lda,B, ldb) 						strsm_(side, trian, trans, diag, &m, &n, &alpha, A, &lda, B, &ldb)

#else

#define BLAS_cp(n, dx, incx, dy, incy) 														dcopy_(&(n), dx, &(incx), dy, &(incy))
#define BLAS_gemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) 			dgemm_(transa, transb, &(m), &(n), &(k), &(alpha), A, &(lda), B, &(ldb), &(beta), C, &(ldc))
#define BLAS_dot(n, dx, incx, dy, incy) 													ddot_(&(n), dx, &(incx), dy, &(incy)) 
#define BLAS_axpy(n, da, dx, incx, dy, incy) 												daxpy_(&n, &da, dx, &incx, dy, &incy)
#define BLAS_scal(n, da, dx, incx)															dscal_(&n, &da, dx, &incx)
#define BLAS_nrm2(n, x, incx)																dnrm2_(&n, x, &incx)
#define BLAS_syrk(trian, trans, n, k, alpha, A, lda, beta, C, ldc)							dsyrk_(trian, trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc)		
#define BLAS_potrf(trian, n, A, lda, info)													dpotrf_(trian, &n, A, &lda, &info)
#define BLAS_trsm(side, trian, trans, diag, m, n, alpha,A, lda,B, ldb) 						dtrsm_(side, trian, trans, diag, &m, &n, &alpha, A, &lda, B, &ldb)

#endif // PRECISION


#endif // MKL


#endif // __FPBLAS_H__
