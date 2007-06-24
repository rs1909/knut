#include "f2c.h"

logical lsame_(char *ca, char *cb, ftnlen ca_len, ftnlen cb_len);

void daxpy_(integer * N, doublereal * alpha, doublereal * X, integer * incX, doublereal * Y, integer * incY);

void dcopy_(integer * N, doublereal * X, integer * incX, doublereal * Y, integer * incY);

doublereal ddot_(integer * N, doublereal * X, integer * incX, doublereal * Y, integer * incY);

void dgemm_(char *TRANSA, char *TRANSB, integer * M, integer * N, integer * K, doublereal * alpha, doublereal * A, integer * lda, doublereal * B, integer * ldb, doublereal * beta, doublereal * C, integer * ldc, ftnlen l1, ftnlen l2);

void dgemv_(char *TRANSA, integer * M, integer * N, doublereal * alpha, doublereal * A, integer * lda, doublereal * X, integer * incX, doublereal * beta, doublereal * Y, integer * incY, ftnlen l1);

void dtrsm_(char *SIDE, char *UPLO, char *TRANSA, char *DIAG, integer * M, integer * N, doublereal * alpha, doublereal * A, integer * lda, doublereal * B, integer * ldb, ftnlen l1, ftnlen l2, ftnlen l3, ftnlen l4);

void dtrmv_(char *UPLO, char *TRANSA, char *DIAG, integer * N, doublereal * A, integer * lda, doublereal * X, integer * incX, ftnlen l1, ftnlen l2, ftnlen l3);

void dtrsv_(char *UPLO, char *TRANSA, char *DIAG, integer * N, doublereal * A, integer * lda, doublereal * X, integer * incX, ftnlen l1, ftnlen l2, ftnlen l3);

void dtrmm_(char *SIDE, char *UPLO, char *TRANSA, char *DIAG, integer * M, integer * N, doublereal * alpha, doublereal * A, integer * lda, doublereal * B, integer * ldb, ftnlen l1, ftnlen l2, ftnlen l3, ftnlen l4);

void dscal_(integer * N, doublereal * alpha, doublereal * X, integer * incX);

void dswap_(integer * N, doublereal * X, integer * incX, doublereal * Y, integer * incY);

integer idamax_(integer * N, doublereal * X, integer * incX);

doublereal dnrm2_(integer * N, doublereal * X, integer * incX);

void drot_(integer * N, doublereal * X, integer * incX, doublereal * Y, integer * incY, doublereal * c, doublereal * s);

doublereal dasum_(integer * N, doublereal * X, integer * incX);

void dger_(integer * M, integer * N, doublereal * alpha, doublereal * X, integer * incX, doublereal * Y, integer * incY, doublereal * A, integer * lda);
