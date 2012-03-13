#ifndef __BLAS_H
#define __BLAS_H

void daxpy_(int * N, double * alpha, double * X, int * incX, double * Y, int * incY);
void dcopy_(int * N, double * X, int * incX, double * Y, int * incY);
double ddot_(int * N, double * X, int * incX, double * Y, int * incY);
void dgemm_(char *TRANSA, char *TRANSB, int * M, int * N, int * K, double * alpha, double * A, int * lda, double * B, int * ldb,
       double * beta, double * C, int * ldc, int l1, int l2);
void dgemv_(char *TRANSA, int * M, int * N, double * alpha, double * A, int * lda, double * X, int * incX, double * beta,
       double * Y, int * incY, int l1);
void dtrsm_(char *SIDE, char *UPLO, char *TRANSA, char *DIAG, int * M, int * N, double * alpha, double * A, int * lda, double * B,
       int * ldb, int l1, int l2, int l3, int l4);
void dtrmv_(char *UPLO, char *TRANSA, char *DIAG, int * N, double * A, int * lda, double * X, int * incX, int l1, int l2, int l3);
void dtrsv_(char *UPLO, char *TRANSA, char *DIAG, int * N, double * A, int * lda, double * X, int * incX, int l1, int l2, int l3);
void dtrmm_(char *SIDE, char *UPLO, char *TRANSA, char *DIAG, int * M, int * N, double * alpha, double * A, int * lda, double * B,
       int * ldb, int l1, int l2, int l3, int l4);
void dscal_(int * N, double * alpha, double * X, int * incX);
void dswap_(int * N, double * X, int * incX, double * Y, int * incY);
int idamax_(int * N, double * X, int * incX);
double dnrm2_(int * N, double * X, int * incX);
void drot_(int * N, double * X, int * incX, double * Y, int * incY, double * c, double * s);
double dasum_(int * N, double * X, int * incX);
void dger_(int * M, int * N, double * alpha, double * X, int * incX, double * Y, int * incY, double * A, int * lda);

// These used in matrix.h
inline void   BLAS_daxpy(size_t N, double alpha, double * X, size_t incX, double * Y, size_t incY)
{
  int N_ = (int)N;
  int incX_ = (int)incX;
  int incY_ = (int)incY;
  daxpy_(&N_, &alpha, X, &incX_, Y, &incY_);
}
inline void   BLAS_dcopy(size_t N, double * X, size_t incX, double * Y, size_t incY)
{
  int N_ = (int)N;
  int incX_ = (int)incX;
  int incY_ = (int)incY;
  dcopy_(&N_, X, &incX_, Y, &incY_);
}
inline double BLAS_ddot(size_t N, double * X, size_t incX, double * Y, size_t incY)
{
  int N_ = (int)N;
  int incX_ = (int)incX;
  int incY_ = (int)incY;
  return ddot_(&N_, X, &incX_, Y, &incY_);
}
inline void   BLAS_dgemm(char TRANSA, char TRANSB, size_t M, size_t N, size_t K, double alpha, double * A, size_t lda, double * B, size_t ldb,
       double beta, double * C, size_t ldc)
{
  int M_ = (int)M;
  int N_ = (int)N;
  int K_ = (int)K;
  int lda_ = (int)lda;
  int ldb_ = (int)ldb;
  int ldc_ = (int)ldc;
  dgemm_(&TRANSA, &TRANSB, &M_, &N_, &K_, &alpha, A, &lda_, B, &ldb_, &beta, C, &ldc_, 1, 1);
}
inline void   BLAS_dgemv(char TRANSA, size_t M, size_t N, double alpha, double * A, size_t lda, double * X, size_t incX, double beta,
       double * Y, size_t incY)
{
  int M_ = (int)M;
  int N_ = (int)N;
  int lda_ = (int)lda;
  int incX_ = (int)incX;
  int incY_ = (int)incY;
  dgemv_(&TRANSA, &M_, &N_, &alpha, A, &lda_, X, &incX_, &beta, Y, &incY_, 1);
}
inline void   BLAS_dscal(size_t N, double alpha, double * X, size_t incX)
{
  int N_ = (int)N;
  int incX_ = (int)incX;
  dscal_(&N_, &alpha, X, &incX_);
}

#endif
