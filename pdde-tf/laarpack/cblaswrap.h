#include "cblas.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

// for defining bool in a c++ compatible way
#include <stdbool.h>

#define max(a,b) ( (a) > b ? (a) : (b) )
#define min(a,b) ( (a) > b ? (b) : (a) )
#define abs(a) (((a)>0)?(a):(-a))

#define s_cmp(a,b,c,d) strncmp(a,b,min(c,d))
#define s_copy(a,b,c,d) strncpy(a,b,min(c,d))

static inline float etime_(float* et)
{
  clock_t clicks;
  clicks = clock();
  return ((float) clicks / CLOCKS_PER_SEC);
}

static inline double d_lg10(double* a)
{
  return log10(*a);
}

static inline double pow_di(double* a, int *i)
{
  return pow(*a, *i);
}

static inline double pow_dd(double* a, double *b)
{
  return pow(*a, *b);
}

static inline double d_sign(double * a, double * b)
{
  double x;
  x = (*a >= 0 ? *a : -*a);
  return (*b >= 0 ? x : -x);
}

static inline int i_dnnt(double * x)
{
  return (int)(*x >= 0. ? floor(*x + .5) : -floor(.5 - *x));
}

static inline char *f77_alloc(int Len, char *whence)
{
  char *rv;
  if (!(rv = (char *) malloc(Len)))
  {
    fprintf(stderr, "malloc(%d) failure in %s\n", Len, whence);
    exit(3);
  }
  return rv;
}

static inline void s_cat(char *lp, char *rpp[], int rnp[], int * np, int ll)
{
  int i, nc;
  char *rp;
  int n = *np;
  int L, m;
  char *lp0, *lp1;

  lp0 = 0;
  lp1 = lp;
  L = ll;
  i = 0;
  while (i < n)
  {
    rp = rpp[i];
    m = rnp[i++];
    if (rp >= lp1 || rp + m <= lp)
    {
      if ((L -= m) <= 0)
      {
        n = i;
        break;
      }
      lp1 += m;
      continue;
    }
    lp0 = lp;
    lp = lp1 = f77_alloc(L = ll, "s_cat");
    break;
  }
  lp1 = lp;
  for (i = 0; i < n; ++i)
  {
    nc = ll;
    if (rnp[i] < nc)
      nc = rnp[i];
    ll -= nc;
    rp = rpp[i];
    while (--nc >= 0)
      *lp++ = *rp++;
  }
  while (--ll >= 0)
    *lp++ = ' ';
  if (lp0)
  {
    memcpy(lp0, lp1, L);
    free(lp1);
  }
}

static inline bool
lsame_(char *ca, char *cb, int ca_len, int cb_len)
{
  if (*ca == *cb)
    return true;
  else
    return toupper(*ca) == toupper(*cb);
}

static inline void
daxpy_(int * N, double * alpha, double * X, int * incX, double * Y, int * incY)
{
  cblas_daxpy(*N, *alpha, X, *incX, Y, *incY);
}

static inline void
dcopy_(int * N, double * X, int * incX, double * Y, int * incY)
{
  cblas_dcopy(*N, X, *incX, Y, *incY);
}

static inline double
ddot_(int * N, double * X, int * incX, double * Y, int * incY)
{
  return cblas_ddot(*N, X, *incX, Y, *incY);
}

static inline void
dgemm_(char *TRANSA, char *TRANSB, int * M, int * N, int * K, double * alpha, double * A, int * lda, double * B, int * ldb,
       double * beta, double * C, int * ldc, int l1, int l2)
{
  enum CBLAS_TRANSPOSE TransA;
  enum CBLAS_TRANSPOSE TransB;

  if (lsame_(TRANSA, "N", 1, 1))
    TransA = CblasNoTrans;
  else if (lsame_(TRANSA, "T", 1, 1))
    TransA = CblasTrans;
  else if (lsame_(TRANSA, "C", 1, 1))
    TransA = CblasConjTrans;
  else
  {
    printf("Error\n");
    return;
  }

  if (lsame_(TRANSB, "N", 1, 1))
    TransB = CblasNoTrans;
  else if (lsame_(TRANSB, "T", 1, 1))
    TransB = CblasTrans;
  else if (lsame_(TRANSB, "C", 1, 1))
    TransB = CblasConjTrans;
  else
  {
    printf("Error\n");
    return;
  }

  cblas_dgemm(CblasColMajor, TransA, TransB, *M, *N, *K, *alpha, A, *lda, B, *ldb, *beta, C, *ldc);
}

static inline void
dgemv_(char *TRANSA, int * M, int * N, double * alpha, double * A, int * lda, double * X, int * incX, double * beta,
       double * Y, int * incY, int l1)
{
  enum CBLAS_TRANSPOSE TransA;

  if (lsame_(TRANSA, "N", 1, 1))
    TransA = CblasNoTrans;
  else if (lsame_(TRANSA, "T", 1, 1))
    TransA = CblasTrans;
  else if (lsame_(TRANSA, "C", 1, 1))
    TransA = CblasConjTrans;
  else
  {
    printf("Error\n");
    return;
  }

  cblas_dgemv(CblasColMajor, TransA, *M, *N, *alpha, A, *lda, X, *incX, *beta, Y, *incY);
}

static inline void
dtrsm_(char *SIDE, char *UPLO, char *TRANSA, char *DIAG, int * M, int * N, double * alpha, double * A, int * lda, double * B,
       int * ldb, int l1, int l2, int l3, int l4)
{
  enum CBLAS_SIDE Side;
  enum CBLAS_UPLO Uplo;
  enum CBLAS_TRANSPOSE TransA;
  enum CBLAS_DIAG Diag;

  if (lsame_(SIDE, "L", 1, 1))
    Side = CblasLeft;
  else if (lsame_(SIDE, "R", 1, 1))
    Side = CblasRight;
  else
  {
    printf("Error\n");
    return;
  }

  if (lsame_(UPLO, "U", 1, 1))
    Uplo = CblasUpper;
  else if (lsame_(UPLO, "L", 1, 1))
    Uplo = CblasLower;
  else
  {
    printf("Error\n");
    return;
  }

  if (lsame_(TRANSA, "N", 1, 1))
    TransA = CblasNoTrans;
  else if (lsame_(TRANSA, "T", 1, 1))
    TransA = CblasTrans;
  else if (lsame_(TRANSA, "C", 1, 1))
    TransA = CblasConjTrans;
  else
  {
    printf("Error\n");
    return;
  }

  if (lsame_(DIAG, "N", 1, 1))
    Diag = CblasNonUnit;
  else if (lsame_(DIAG, "U", 1, 1))
    Diag = CblasUnit;
  else
  {
    printf("Error\n");
    return;
  }

  cblas_dtrsm(CblasColMajor, Side, Uplo, TransA, Diag, *M, *N, *alpha, A, *lda, B, *ldb);
}

static inline void
dtrmv_(char *UPLO, char *TRANSA, char *DIAG, int * N, double * A, int * lda, double * X, int * incX, int l1, int l2, int l3)
{
  enum CBLAS_UPLO Uplo;
  enum CBLAS_TRANSPOSE TransA;
  enum CBLAS_DIAG Diag;

  if (lsame_(UPLO, "U", 1, 1))
    Uplo = CblasUpper;
  else if (lsame_(UPLO, "L", 1, 1))
    Uplo = CblasLower;
  else
  {
    printf("Error\n");
    return;
  }

  if (lsame_(TRANSA, "N", 1, 1))
    TransA = CblasNoTrans;
  else if (lsame_(TRANSA, "T", 1, 1))
    TransA = CblasTrans;
  else if (lsame_(TRANSA, "C", 1, 1))
    TransA = CblasConjTrans;
  else
  {
    printf("Error\n");
    return;
  }

  if (lsame_(DIAG, "N", 1, 1))
    Diag = CblasNonUnit;
  else if (lsame_(DIAG, "U", 1, 1))
    Diag = CblasUnit;
  else
  {
    printf("Error\n");
    return;
  }

  cblas_dtrmv(CblasColMajor, Uplo, TransA, Diag, *N, A, *lda, X, *incX);
}

static inline void
dtrsv_(char *UPLO, char *TRANSA, char *DIAG, int * N, double * A, int * lda, double * X, int * incX, int l1, int l2, int l3)
{
  enum CBLAS_UPLO Uplo;
  enum CBLAS_TRANSPOSE TransA;
  enum CBLAS_DIAG Diag;

  if (lsame_(UPLO, "U", 1, 1))
    Uplo = CblasUpper;
  else if (lsame_(UPLO, "L", 1, 1))
    Uplo = CblasLower;
  else
  {
    printf("Error\n");
    return;
  }

  if (lsame_(TRANSA, "N", 1, 1))
    TransA = CblasNoTrans;
  else if (lsame_(TRANSA, "T", 1, 1))
    TransA = CblasTrans;
  else if (lsame_(TRANSA, "C", 1, 1))
    TransA = CblasConjTrans;
  else
  {
    printf("Error\n");
    return;
  }

  if (lsame_(DIAG, "N", 1, 1))
    Diag = CblasNonUnit;
  else if (lsame_(DIAG, "U", 1, 1))
    Diag = CblasUnit;
  else
  {
    printf("Error\n");
    return;
  }

  cblas_dtrsv(CblasColMajor, Uplo, TransA, Diag, *N, A, *lda, X, *incX);
}

static inline void
dtrmm_(char *SIDE, char *UPLO, char *TRANSA, char *DIAG, int * M, int * N, double * alpha, double * A, int * lda, double * B,
       int * ldb, int l1, int l2, int l3, int l4)
{
  enum CBLAS_SIDE Side;
  enum CBLAS_UPLO Uplo;
  enum CBLAS_TRANSPOSE TransA;
  enum CBLAS_DIAG Diag;

  if (lsame_(SIDE, "L", 1, 1))
    Side = CblasLeft;
  else if (lsame_(SIDE, "R", 1, 1))
    Side = CblasRight;
  else
  {
    printf("Error\n");
    return;
  }

  if (lsame_(UPLO, "U", 1, 1))
    Uplo = CblasUpper;
  else if (lsame_(UPLO, "L", 1, 1))
    Uplo = CblasLower;
  else
  {
    printf("Error\n");
    return;
  }

  if (lsame_(TRANSA, "N", 1, 1))
    TransA = CblasNoTrans;
  else if (lsame_(TRANSA, "T", 1, 1))
    TransA = CblasTrans;
  else if (lsame_(TRANSA, "C", 1, 1))
    TransA = CblasConjTrans;
  else
  {
    printf("Error\n");
    return;
  }

  if (lsame_(DIAG, "N", 1, 1))
    Diag = CblasNonUnit;
  else if (lsame_(DIAG, "U", 1, 1))
    Diag = CblasUnit;
  else
  {
    printf("Error\n");
    return;
  }

  cblas_dtrmm(CblasColMajor, Side, Uplo, TransA, Diag, *M, *N, *alpha, A, *lda, B, *ldb);
}

static inline void
dscal_(int * N, double * alpha, double * X, int * incX)
{
  cblas_dscal(*N, *alpha, X, *incX);
}

static inline void
dswap_(int * N, double * X, int * incX, double * Y, int * incY)
{
  cblas_dswap(*N, X, *incX, Y, *incY);
}

static inline int
idamax_(int * N, double * X, int * incX)
{
  return cblas_idamax(*N, X, *incX) + 1;
}

static inline double
dnrm2_(int * N, double * X, int * incX)
{
  return cblas_dnrm2(*N, X, *incX);
}

static inline void
drot_(int * N, double * X, int * incX, double * Y, int * incY, double * c, double * s)
{
  cblas_drot(*N, X, *incX, Y, *incY, *c, *s);
}

static inline double
dasum_(int * N, double * X, int * incX)
{
  return cblas_dasum(*N, X, *incX);
}

static inline void
dger_(int * M, int * N, double * alpha, double * X, int * incX, double * Y, int * incY, double * A, int * lda)
{
  cblas_dger(CblasColMajor, *M, *N, *alpha, X, *incX, Y, *incY, A, *lda);
}
