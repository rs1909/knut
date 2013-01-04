#ifndef __BLAS_H
#define __BLAS_H

#include "config.h"
#include <stddef.h>

#ifdef LONGBLAS
typedef ptrdiff_t blasint;
#else
typedef int blasint;
#endif

void daxpy_(blasint * N, double * alpha, double * X, blasint * incX, double * Y, blasint * incY);
void dcopy_(blasint * N, double * X, blasint * incX, double * Y, blasint * incY);
double ddot_(blasint * N, double * X, blasint * incX, double * Y, blasint * incY);
void dgemm_(char *TRANSA, char *TRANSB, blasint * M, blasint * N, blasint * K, double * alpha, double * A, blasint * lda, double * B, blasint * ldb,
       double * beta, double * C, blasint * ldc, blasint l1, blasint l2);
void dgemv_(char *TRANSA, blasint * M, blasint * N, double * alpha, double * A, blasint * lda, double * X, blasint * incX, double * beta,
       double * Y, blasint * incY, blasint l1);
void dtrsm_(char *SIDE, char *UPLO, char *TRANSA, char *DIAG, blasint * M, blasint * N, double * alpha, double * A, blasint * lda, double * B,
       blasint * ldb, blasint l1, blasint l2, blasint l3, blasint l4);
void dtrmv_(char *UPLO, char *TRANSA, char *DIAG, blasint * N, double * A, blasint * lda, double * X, blasint * incX, blasint l1, blasint l2, blasint l3);
void dtrsv_(char *UPLO, char *TRANSA, char *DIAG, blasint * N, double * A, blasint * lda, double * X, blasint * incX, blasint l1, blasint l2, blasint l3);
void dtrmm_(char *SIDE, char *UPLO, char *TRANSA, char *DIAG, blasint * M, blasint * N, double * alpha, double * A, blasint * lda, double * B,
       blasint * ldb, blasint l1, blasint l2, blasint l3, blasint l4);
void dscal_(blasint * N, double * alpha, double * X, blasint * incX);
void dswap_(blasint * N, double * X, blasint * incX, double * Y, blasint * incY);
blasint idamax_(blasint * N, double * X, blasint * incX);
double dnrm2_(blasint * N, double * X, blasint * incX);
void drot_(blasint * N, double * X, blasint * incX, double * Y, blasint * incY, double * c, double * s);
double dasum_(blasint * N, double * X, blasint * incX);
void dger_(blasint * M, blasint * N, double * alpha, double * X, blasint * incX, double * Y, blasint * incY, double * A, blasint * lda);

// These used in matrix.h
inline void   BLAS_daxpy(ptrdiff_t N, double alpha, double * X, ptrdiff_t incX, double * Y, ptrdiff_t incY)
{
  blasint N_ = (blasint)N;
  blasint incX_ = (blasint)incX;
  blasint incY_ = (blasint)incY;
  daxpy_(&N_, &alpha, X, &incX_, Y, &incY_);
}
inline void   BLAS_dcopy(ptrdiff_t N, double * X, ptrdiff_t incX, double * Y, ptrdiff_t incY)
{
  blasint N_ = (blasint)N;
  blasint incX_ = (blasint)incX;
  blasint incY_ = (blasint)incY;
  dcopy_(&N_, X, &incX_, Y, &incY_);
}
inline double BLAS_ddot(ptrdiff_t N, double * X, ptrdiff_t incX, double * Y, ptrdiff_t incY)
{
  blasint N_ = (blasint)N;
  blasint incX_ = (blasint)incX;
  blasint incY_ = (blasint)incY;
  return ddot_(&N_, X, &incX_, Y, &incY_);
}
inline void   BLAS_dgemm(char TRANSA, char TRANSB, ptrdiff_t M, ptrdiff_t N, ptrdiff_t K, double alpha, double * A, ptrdiff_t lda, double * B, ptrdiff_t ldb,
       double beta, double * C, ptrdiff_t ldc)
{
  blasint M_ = (blasint)M;
  blasint N_ = (blasint)N;
  blasint K_ = (blasint)K;
  blasint lda_ = (blasint)lda;
  blasint ldb_ = (blasint)ldb;
  blasint ldc_ = (blasint)ldc;
  dgemm_(&TRANSA, &TRANSB, &M_, &N_, &K_, &alpha, A, &lda_, B, &ldb_, &beta, C, &ldc_, 1, 1);
}
inline void   BLAS_dgemv(char TRANSA, ptrdiff_t M, ptrdiff_t N, double alpha, double * A, ptrdiff_t lda, double * X, ptrdiff_t incX, double beta,
       double * Y, ptrdiff_t incY)
{
  blasint M_ = (blasint)M;
  blasint N_ = (blasint)N;
  blasint lda_ = (blasint)lda;
  blasint incX_ = (blasint)incX;
  blasint incY_ = (blasint)incY;
  dgemv_(&TRANSA, &M_, &N_, &alpha, A, &lda_, X, &incX_, &beta, Y, &incY_, 1);
}
inline void   BLAS_dscal(ptrdiff_t N, double alpha, double * X, ptrdiff_t incX)
{
  blasint N_ = (blasint)N;
  blasint incX_ = (blasint)incX;
  dscal_(&N_, &alpha, X, &incX_);
}

#endif
