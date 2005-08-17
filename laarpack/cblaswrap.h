#include "f2c.h"
#include "cblas.h"
#include <ctype.h>
#include <stdio.h>

static inline logical lsame_(char *ca,char *cb,ftnlen ca_len,ftnlen cb_len)
{
	if( *ca == *cb ) return TRUE_;
	else return toupper( *ca ) == toupper( *cb );
}

static inline int daxpy_(integer* N, doublereal* alpha, doublereal* X, integer* incX, doublereal* Y, integer* incY)
{
	cblas_daxpy( *N, *alpha, X, *incX, Y, *incY );
}

static inline int  dcopy_(integer* N, doublereal* X, integer* incX, doublereal* Y, integer* incY)
{
	cblas_dcopy( *N, X, *incX, Y, *incY );
}

static inline doublereal ddot_(integer* N, doublereal* X, integer* incX, doublereal* Y, integer* incY)
{
	return cblas_ddot( *N, X, *incX, Y, *incY );
}

static inline int  dgemm_(char *TRANSA, char *TRANSB, integer *M, integer *N, integer *K, doublereal *alpha, doublereal *A, integer *lda, doublereal *B, integer *ldb, doublereal *beta, doublereal *C, integer *ldc, ftnlen l1, ftnlen l2)
{
	enum CBLAS_TRANSPOSE TransA;
	enum CBLAS_TRANSPOSE TransB;

	if( lsame_( TRANSA, "N", 1, 1 ) ) TransA = CblasNoTrans;
	else if( lsame_( TRANSA, "T", 1, 1 ) ) TransA = CblasTrans;
	else if( lsame_( TRANSA, "C", 1, 1 ) ) TransA = CblasConjTrans;
	else printf("Error\n");

	if( lsame_( TRANSB, "N", 1, 1 ) ) TransB = CblasNoTrans;
	else if( lsame_( TRANSB, "T", 1, 1 ) ) TransB = CblasTrans;
	else if( lsame_( TRANSB, "C", 1, 1 ) ) TransB = CblasConjTrans;
	else printf("Error\n");

	cblas_dgemm( CblasColMajor, TransA, TransB, *M,*N, *K, *alpha, A, *lda, B, *ldb, *beta, C, *ldc);
}

static inline int  dgemv_(char *TRANSA, integer* M, integer* N, doublereal* alpha, doublereal* A, integer* lda, doublereal* X, integer* incX, doublereal* beta, doublereal* Y, integer* incY, ftnlen l1)
{
	enum CBLAS_TRANSPOSE TransA;

	if( lsame_( TRANSA, "N", 1, 1 ) ) TransA = CblasNoTrans;
	else if( lsame_( TRANSA, "T", 1, 1 ) ) TransA = CblasTrans;
	else if( lsame_( TRANSA, "C", 1, 1 ) ) TransA = CblasConjTrans;
	else printf("Error\n");

	cblas_dgemv( CblasColMajor, TransA, *M, *N, *alpha, A, *lda, X, *incX, *beta, Y, *incY);
}

static inline int dtrsm_(char *SIDE, char *UPLO, char *TRANSA, char *DIAG, integer *M, integer *N, doublereal *alpha, doublereal *A, integer *lda, doublereal *B, integer *ldb, ftnlen l1, ftnlen l2, ftnlen l3, ftnlen l4)
{
	enum CBLAS_SIDE Side;
	enum CBLAS_UPLO Uplo;
	enum CBLAS_TRANSPOSE TransA;
	enum CBLAS_DIAG Diag;

	if( lsame_( SIDE, "L", 1, 1 ) ) Side = CblasLeft;
	else if( lsame_( SIDE, "R", 1, 1 ) ) Side = CblasRight;
	else printf("Error\n");

	if( lsame_( UPLO, "U", 1, 1 ) ) Uplo = CblasUpper;
	else if( lsame_( UPLO, "L", 1, 1 ) ) Uplo = CblasLower;
	else printf("Error\n");

	if( lsame_( TRANSA, "N", 1, 1 ) ) TransA = CblasNoTrans;
	else if( lsame_( TRANSA, "T", 1, 1 ) ) TransA = CblasTrans;
	else if( lsame_( TRANSA, "C", 1, 1 ) ) TransA = CblasConjTrans;
	else printf("Error\n");

	if( lsame_( DIAG, "N", 1, 1 ) ) Diag = CblasNonUnit;
	else if( lsame_( DIAG, "U", 1, 1 ) ) Diag = CblasUnit;
	else printf("Error\n");

	cblas_dtrsm( CblasColMajor, Side, Uplo, TransA, Diag, *M, *N, *alpha, A, *lda, B, *ldb);
}

static inline int dtrmv_(char *UPLO, char *TRANSA, char *DIAG, integer *N, doublereal *A, integer *lda, doublereal *X, integer *incX, ftnlen l1, ftnlen l2, ftnlen l3)
{
	enum CBLAS_UPLO Uplo;
	enum CBLAS_TRANSPOSE TransA;
	enum CBLAS_DIAG Diag;

	if( lsame_( UPLO, "U", 1, 1 ) ) Uplo = CblasUpper;
	else if( lsame_( UPLO, "L", 1, 1 ) ) Uplo = CblasLower;
	else printf("Error\n");

	if( lsame_( TRANSA, "N", 1, 1 ) ) TransA = CblasNoTrans;
	else if( lsame_( TRANSA, "T", 1, 1 ) ) TransA = CblasTrans;
	else if( lsame_( TRANSA, "C", 1, 1 ) ) TransA = CblasConjTrans;
	else printf("Error\n");

	if( lsame_( DIAG, "N", 1, 1 ) ) Diag = CblasNonUnit;
	else if( lsame_( DIAG, "U", 1, 1 ) ) Diag = CblasUnit;
	else printf("Error\n");

	cblas_dtrmv( CblasColMajor, Uplo, TransA, Diag, *N, A, *lda, X, *incX );
}

static inline int dtrsv_(char *UPLO, char *TRANSA, char *DIAG, integer *N, doublereal *A, integer *lda, doublereal *X, integer *incX, ftnlen l1, ftnlen l2, ftnlen l3)
{
	enum CBLAS_UPLO Uplo;
	enum CBLAS_TRANSPOSE TransA;
	enum CBLAS_DIAG Diag;

	if( lsame_( UPLO, "U", 1, 1 ) ) Uplo = CblasUpper;
	else if( lsame_( UPLO, "L", 1, 1 ) ) Uplo = CblasLower;
	else printf("Error\n");

	if( lsame_( TRANSA, "N", 1, 1 ) ) TransA = CblasNoTrans;
	else if( lsame_( TRANSA, "T", 1, 1 ) ) TransA = CblasTrans;
	else if( lsame_( TRANSA, "C", 1, 1 ) ) TransA = CblasConjTrans;
	else printf("Error\n");

	if( lsame_( DIAG, "N", 1, 1 ) ) Diag = CblasNonUnit;
	else if( lsame_( DIAG, "U", 1, 1 ) ) Diag = CblasUnit;
	else printf("Error\n");

	cblas_dtrsv( CblasColMajor, Uplo, TransA, Diag, *N, A, *lda, X, *incX );
}

static inline int dtrmm_(char *SIDE, char *UPLO, char *TRANSA, char *DIAG, integer *M, integer *N, doublereal *alpha, doublereal *A, integer *lda, doublereal *B, integer *ldb, ftnlen l1, ftnlen l2, ftnlen l3, ftnlen l4)
{
	enum CBLAS_SIDE Side;
	enum CBLAS_UPLO Uplo;
	enum CBLAS_TRANSPOSE TransA;
	enum CBLAS_DIAG Diag;

	if( lsame_( SIDE, "L", 1, 1 ) ) Side = CblasLeft;
	else if( lsame_( SIDE, "R", 1, 1 ) ) Side = CblasRight;
	else printf("Error\n");

	if( lsame_( UPLO, "U", 1, 1 ) ) Uplo = CblasUpper;
	else if( lsame_( UPLO, "L", 1, 1 ) ) Uplo = CblasLower;
	else printf("Error\n");

	if( lsame_( TRANSA, "N", 1, 1 ) ) TransA = CblasNoTrans;
	else if( lsame_( TRANSA, "T", 1, 1 ) ) TransA = CblasTrans;
	else if( lsame_( TRANSA, "C", 1, 1 ) ) TransA = CblasConjTrans;
	else printf("Error\n");

	if( lsame_( DIAG, "N", 1, 1 ) ) Diag = CblasNonUnit;
	else if( lsame_( DIAG, "U", 1, 1 ) ) Diag = CblasUnit;
	else printf("Error\n");

	cblas_dtrmm( CblasColMajor, Side, Uplo, TransA, Diag, *M, *N, *alpha, A, *lda, B, *ldb );
}

static inline int dscal_(integer* N, doublereal* alpha, doublereal* X, integer* incX)
{
	cblas_dscal( *N, *alpha, X, *incX );
}

static inline int dswap_(integer *N, doublereal *X, integer *incX, doublereal *Y, integer *incY)
{
	cblas_dswap( *N, X, *incX, Y, *incY );
}

static inline integer idamax_(integer *N, doublereal *X, integer *incX)
{
	return cblas_idamax( *N, X, *incX ) + 1;
}

static inline doublereal dnrm2_(integer *N, doublereal *X, integer *incX)
{
	return cblas_dnrm2( *N, X, *incX );
}

static inline int drot_(integer *N, doublereal *X, integer *incX, doublereal *Y, integer *incY, doublereal *c, doublereal *s)
{
	cblas_drot( *N, X, *incX, Y, *incY, *c, *s );
}

static inline doublereal dasum_(integer *N, doublereal *X, integer *incX)
{
	return cblas_dasum( *N, X, *incX );
}

static inline int dger_( integer *M, integer *N, doublereal *alpha, doublereal *X, integer *incX, doublereal *Y, integer *incY, doublereal *A, integer *lda )
{
	cblas_dger( CblasColMajor, *M, *N, *alpha, X, *incX, Y, *incY, A, *lda);
}
