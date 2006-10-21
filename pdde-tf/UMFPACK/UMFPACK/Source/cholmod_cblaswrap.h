#include "cblas.h"
#include <ctype.h>
#include <stdio.h>

typedef int bool;

#define true (1)
#define false (0)

static inline bool
lsame_ (char *ca, char *cb)
{
	if (*ca == *cb)
		return true;
	else
		return toupper (*ca) == toupper (*cb);
}

static inline void
daxpy_ (BLAS_INT * N, double * alpha, double * X, BLAS_INT * incX, double * Y, BLAS_INT * incY)
{
	cblas_daxpy (*N, *alpha, X, *incX, Y, *incY);
}

static inline void
dcopy_ (BLAS_INT * N, double * X, BLAS_INT * incX, double * Y, BLAS_INT * incY)
{
	cblas_dcopy (*N, X, *incX, Y, *incY);
}

static inline double
ddot_ (BLAS_INT * N, double * X, BLAS_INT * incX, double * Y, BLAS_INT * incY)
{
	return cblas_ddot (*N, X, *incX, Y, *incY);
}

static inline void
dgemm_ (char *TRANSA, char *TRANSB, BLAS_INT * M, BLAS_INT * N, BLAS_INT * K, double * alpha, double * A, BLAS_INT * lda, double * B, BLAS_INT * ldb,
		  double * beta, double * C, BLAS_INT * ldc)
{
	enum CBLAS_TRANSPOSE TransA;
	enum CBLAS_TRANSPOSE TransB;

	if (lsame_ (TRANSA, "N"))
		TransA = CblasNoTrans;
	else if (lsame_ (TRANSA, "T"))
		TransA = CblasTrans;
	else if (lsame_ (TRANSA, "C"))
		TransA = CblasConjTrans;
	else
		{ printf ("Error\n"); return; }

	if (lsame_ (TRANSB, "N"))
		TransB = CblasNoTrans;
	else if (lsame_ (TRANSB, "T"))
		TransB = CblasTrans;
	else if (lsame_ (TRANSB, "C"))
		TransB = CblasConjTrans;
	else
		{ printf ("Error\n"); return; }

	cblas_dgemm (CblasColMajor, TransA, TransB, *M, *N, *K, *alpha, A, *lda, B, *ldb, *beta, C, *ldc);
}

static inline void
dgemv_ (char *TRANSA, BLAS_INT * M, BLAS_INT * N, double * alpha, double * A, BLAS_INT * lda, double * X, BLAS_INT * incX, double * beta,
		  double * Y, BLAS_INT * incY)
{
	enum CBLAS_TRANSPOSE TransA;

	if (lsame_ (TRANSA, "N"))
		TransA = CblasNoTrans;
	else if (lsame_ (TRANSA, "T"))
		TransA = CblasTrans;
	else if (lsame_ (TRANSA, "C"))
		TransA = CblasConjTrans;
	else
		{ printf ("Error\n"); return; }

	cblas_dgemv (CblasColMajor, TransA, *M, *N, *alpha, A, *lda, X, *incX, *beta, Y, *incY);
}

static inline void
dtrsm_ (char *SIDE, char *UPLO, char *TRANSA, char *DIAG, BLAS_INT * M, BLAS_INT * N, double * alpha, double * A, BLAS_INT * lda, double * B,
		  BLAS_INT * ldb)
{
	enum CBLAS_SIDE Side;
	enum CBLAS_UPLO Uplo;
	enum CBLAS_TRANSPOSE TransA;
	enum CBLAS_DIAG Diag;

	if (lsame_ (SIDE, "L"))
		Side = CblasLeft;
	else if (lsame_ (SIDE, "R"))
		Side = CblasRight;
	else
		{ printf ("Error\n"); return; }

	if (lsame_ (UPLO, "U"))
		Uplo = CblasUpper;
	else if (lsame_ (UPLO, "L"))
		Uplo = CblasLower;
	else
		{ printf ("Error\n"); return; }

	if (lsame_ (TRANSA, "N"))
		TransA = CblasNoTrans;
	else if (lsame_ (TRANSA, "T"))
		TransA = CblasTrans;
	else if (lsame_ (TRANSA, "C"))
		TransA = CblasConjTrans;
	else
		{ printf ("Error\n"); return; }

	if (lsame_ (DIAG, "N"))
		Diag = CblasNonUnit;
	else if (lsame_ (DIAG, "U"))
		Diag = CblasUnit;
	else
		{ printf ("Error\n"); return; }

	cblas_dtrsm (CblasColMajor, Side, Uplo, TransA, Diag, *M, *N, *alpha, A, *lda, B, *ldb);
}

static inline void
dtrmv_ (char *UPLO, char *TRANSA, char *DIAG, BLAS_INT * N, double * A, BLAS_INT * lda, double * X, BLAS_INT * incX)
{
	enum CBLAS_UPLO Uplo;
	enum CBLAS_TRANSPOSE TransA;
	enum CBLAS_DIAG Diag;

	if (lsame_ (UPLO, "U"))
		Uplo = CblasUpper;
	else if (lsame_ (UPLO, "L"))
		Uplo = CblasLower;
	else
		{ printf ("Error\n"); return; }

	if (lsame_ (TRANSA, "N"))
		TransA = CblasNoTrans;
	else if (lsame_ (TRANSA, "T"))
		TransA = CblasTrans;
	else if (lsame_ (TRANSA, "C"))
		TransA = CblasConjTrans;
	else
		{ printf ("Error\n"); return; }

	if (lsame_ (DIAG, "N"))
		Diag = CblasNonUnit;
	else if (lsame_ (DIAG, "U"))
		Diag = CblasUnit;
	else
		{ printf ("Error\n"); return; }

	cblas_dtrmv (CblasColMajor, Uplo, TransA, Diag, *N, A, *lda, X, *incX);
}

static inline void
dtrsv_ (char *UPLO, char *TRANSA, char *DIAG, BLAS_INT * N, double * A, BLAS_INT * lda, double * X, BLAS_INT * incX)
{
	enum CBLAS_UPLO Uplo;
	enum CBLAS_TRANSPOSE TransA;
	enum CBLAS_DIAG Diag;

	if (lsame_ (UPLO, "U"))
		Uplo = CblasUpper;
	else if (lsame_ (UPLO, "L"))
		Uplo = CblasLower;
	else
		{ printf ("Error\n"); return; }

	if (lsame_ (TRANSA, "N"))
		TransA = CblasNoTrans;
	else if (lsame_ (TRANSA, "T"))
		TransA = CblasTrans;
	else if (lsame_ (TRANSA, "C"))
		TransA = CblasConjTrans;
	else
		{ printf ("Error\n"); return; }

	if (lsame_ (DIAG, "N"))
		Diag = CblasNonUnit;
	else if (lsame_ (DIAG, "U"))
		Diag = CblasUnit;
	else
		{ printf ("Error\n"); return; }

	cblas_dtrsv (CblasColMajor, Uplo, TransA, Diag, *N, A, *lda, X, *incX);
}

static inline void
dtrmm_ (char *SIDE, char *UPLO, char *TRANSA, char *DIAG, BLAS_INT * M, BLAS_INT * N, double * alpha, double * A, BLAS_INT * lda, double * B,
		  BLAS_INT * ldb)
{
	enum CBLAS_SIDE Side;
	enum CBLAS_UPLO Uplo;
	enum CBLAS_TRANSPOSE TransA;
	enum CBLAS_DIAG Diag;

	if (lsame_ (SIDE, "L"))
		Side = CblasLeft;
	else if (lsame_ (SIDE, "R"))
		Side = CblasRight;
	else
		{ printf ("Error\n"); return; }

	if (lsame_ (UPLO, "U"))
		Uplo = CblasUpper;
	else if (lsame_ (UPLO, "L"))
		Uplo = CblasLower;
	else
		{ printf ("Error\n"); return; }

	if (lsame_ (TRANSA, "N"))
		TransA = CblasNoTrans;
	else if (lsame_ (TRANSA, "T"))
		TransA = CblasTrans;
	else if (lsame_ (TRANSA, "C"))
		TransA = CblasConjTrans;
	else
		{ printf ("Error\n"); return; }

	if (lsame_ (DIAG, "N"))
		Diag = CblasNonUnit;
	else if (lsame_ (DIAG, "U"))
		Diag = CblasUnit;
	else
		{ printf ("Error\n"); return; }

	cblas_dtrmm (CblasColMajor, Side, Uplo, TransA, Diag, *M, *N, *alpha, A, *lda, B, *ldb);
}

static inline void
dscal_ (BLAS_INT * N, double * alpha, double * X, BLAS_INT * incX)
{
	cblas_dscal (*N, *alpha, X, *incX);
}

static inline void
dswap_ (BLAS_INT * N, double * X, BLAS_INT * incX, double * Y, BLAS_INT * incY)
{
	cblas_dswap (*N, X, *incX, Y, *incY);
}

static inline BLAS_INT
idamax_ (BLAS_INT * N, double * X, BLAS_INT * incX)
{
	return cblas_idamax (*N, X, *incX) + 1;
}

static inline double
dnrm2_ (BLAS_INT * N, double * X, BLAS_INT * incX)
{
	return cblas_dnrm2 (*N, X, *incX);
}

static inline void
drot_ (BLAS_INT * N, double * X, BLAS_INT * incX, double * Y, BLAS_INT * incY, double * c, double * s)
{
	cblas_drot (*N, X, *incX, Y, *incY, *c, *s);
}

static inline double
dasum_ (BLAS_INT * N, double * X, BLAS_INT * incX)
{
	return cblas_dasum (*N, X, *incX);
}

static inline void
dger_ (BLAS_INT * M, BLAS_INT * N, double * alpha, double * X, BLAS_INT * incX, double * Y, BLAS_INT * incY, double * A, BLAS_INT * lda)
{
	cblas_dger (CblasColMajor, *M, *N, *alpha, X, *incX, Y, *incY, A, *lda);
}
