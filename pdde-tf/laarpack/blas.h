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
inline void   BLAS_daxpy(int N, double alpha, double * X, int incX, double * Y, int incY)
{
  daxpy_(&N, &alpha, X, &incX, Y, &incY);
}
inline void   BLAS_dcopy(int N, double * X, int incX, double * Y, int incY)
{
  dcopy_(&N, X, &incX, Y, &incY);
}
inline double BLAS_ddot(int N, double * X, int incX, double * Y, int incY)
{
  return ddot_(&N, X, &incX, Y, &incY);
}
inline void   BLAS_dgemm(char TRANSA, char TRANSB, int M, int N, int K, double alpha, double * A, int lda, double * B, int ldb,
       double beta, double * C, int ldc)
{
  dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &lda, B, &ldb, &beta, C, &ldc, 1, 1);
}
inline void   BLAS_dgemv(char TRANSA, int M, int N, double alpha, double * A, int lda, double * X, int incX, double beta,
       double * Y, int incY)
{
  dgemv_(&TRANSA, &M, &N, &alpha, A, &lda, X, &incX, &beta, Y, &incY, 1);
}
inline void   BLAS_dscal(int N, double alpha, double * X, int incX)
{
  dscal_(&N, &alpha, X, &incX);
}

#endif
