#ifndef LAARPACK_H
#define LAARPACK_H

#include "config.h"
#include <stddef.h>

#ifdef LONGBLAS
typedef ptrdiff_t blasint;
#else
typedef blasint blasint;
#endif

blasint dneupd_(blasint * rvec, char *howmny, bool * select, double *dr, double *di,
	     double *z__, blasint *ldz, double *sigmar, double *sigmai,
	     double *workev, char *bmat, blasint *n, char *which, blasint *nev,
	     double *tol, double *resid, blasint *ncv, double *v, blasint *ldv,
	     blasint *iparam, blasint *ipntr, double *workd, double *workl, blasint *lworkl,
	     blasint *info, blasint howmny_len, blasint bmat_len, blasint which_len);

blasint dnaupd_(blasint *ido, char *bmat, blasint *n, char *which, blasint *nev, double *tol,
	     double *resid, blasint *ncv, double *v, blasint *ldv, blasint *iparam,
	     blasint *ipntr, double *workd, double *workl, blasint *lworkl, blasint *info,
	     blasint bmat_len, blasint which_len);

blasint dgesvx_(char *fact, char *trans, blasint *n, blasint *nrhs, double *a, blasint *lda,
	     double *af, blasint *ldaf, blasint *ipiv, char *equed, double *r__,
	     double *c__, double *b, blasint *ldb, double *x, blasint *ldx,
	     double *rcond, double *ferr, double *berr, double *work,
	     blasint *iwork, blasint *info, blasint fact_len, blasint trans_len, blasint equed_len);

blasint dgesv_(blasint *n, blasint *nrhs, double *a, blasint *lda, blasint *ipiv, double *b,
	    blasint *ldb, blasint *info);

blasint dgeevx_(char *balanc, char *jobvl, char *jobvr, char *sense, blasint *n,
	     double *a, blasint *lda, double *wr, double *wi, double *vl, blasint *ldvl,
	     double *vr, blasint *ldvr, blasint *ilo, blasint *ihi, double *scale,
	     double *abnrm, double *rconde, double *rcondv, double *work,
	     blasint *lwork, blasint *iwork, blasint *info, blasint balanc_len, blasint jobvl_len,
	     blasint jobvr_len, blasint sense_len);

blasint dgeev_(char *jobvl, char *jobvr, blasint *n, double *a, blasint *lda, double *wr,
	    double *wi, double *vl, blasint *ldvl, double *vr, blasint *ldvr,
	    double *work, blasint *lwork, blasint *info, blasint jobvl_len, blasint jobvr_len);

blasint dlarnv_(blasint *idist, blasint *iseed, blasint *n, double *x);

double dlamch_(const char *cmach, blasint cmach_len);

// for matrix.cpp and spmatrix.cpp

inline blasint knut_dneupd (bool rvec, char howmny, bool * select, double *dr, double *di,
	     double *z__, ptrdiff_t ldz, double *sigmar, double *sigmai,
	     double *workev, char bmat, ptrdiff_t n, char which1, char which2, ptrdiff_t nev,
	     double *tol, double *resid, ptrdiff_t ncv, double *v, ptrdiff_t ldv,
	     blasint *iparam, blasint *ipntr, double *workd, double *workl, ptrdiff_t lworkl,
	     blasint *info)
{
  blasint n_ = (blasint)n;
  blasint rvec_ = (blasint)rvec;
  blasint ldz_ = (blasint)ldz;
  blasint nev_ = (blasint)nev;  
  blasint ncv_ = (blasint)ncv;  
  blasint ldv_ = (blasint)ldv;
  blasint lworkl_ = (blasint)lworkl;
  char which[2] = {which1, which2};  
  return dneupd_(&rvec_, &howmny, select, dr, di,
    z__, &ldz_, sigmar, sigmai,
    workev, &bmat, &n_, which, &nev_,
    tol, resid, &ncv_, v, &ldv_,
    iparam, ipntr, workd, workl, &lworkl_,
    info, 1, 1, 2);
}

inline blasint knut_dnaupd (blasint *ido, char bmat, ptrdiff_t n, char which1, char which2, ptrdiff_t nev, double *tol,
	     double *resid, ptrdiff_t ncv, double *v, ptrdiff_t ldv, blasint *iparam,
	     blasint *ipntr, double *workd, double *workl, ptrdiff_t lworkl, blasint *info)
{
  blasint n_ = (blasint)n;
  blasint nev_ = (blasint)nev;  
  blasint ncv_ = (blasint)ncv;  
  blasint ldv_ = (blasint)ldv;
  blasint lworkl_ = (blasint)lworkl;
  char which[2] = {which1, which2}; 
  return dnaupd_(ido, &bmat, &n_, which, &nev_, tol,
    resid, &ncv_, v, &ldv_, iparam,
    ipntr, workd, workl, &lworkl_, info,
    1, 2);
}

inline blasint knut_dgesvx (char fact, char trans, ptrdiff_t n, ptrdiff_t nrhs, double *a, ptrdiff_t lda,
	     double *af, ptrdiff_t ldaf, blasint *ipiv, char equed, double *r__,
	     double *c__, double *b, ptrdiff_t ldb, double *x, ptrdiff_t ldx,
	     double *rcond, double *ferr, double *berr, double *work,
	     blasint *iwork, blasint *info)
{
  blasint n_ = (blasint)n;
  blasint nrhs_ = (blasint)nrhs;
  blasint lda_ = (blasint)lda;
  blasint ldb_ = (blasint)ldb;
  blasint ldx_ = (blasint)ldx;
  blasint ldaf_ = (blasint)ldaf;
  return dgesvx_(&fact, &trans, &n_, &nrhs_, a, &lda_,
    af, &ldaf_, ipiv, &equed, r__,
    c__, b, &ldb_, x, &ldx_,
    rcond, ferr, berr, work,
    iwork, info, 1, 1, 1);
}

inline blasint knut_dgesv (blasint *n, blasint *nrhs, double *a, blasint *lda, blasint *ipiv, double *b,
	    blasint *ldb, blasint *info)
{
  return dgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
}

inline blasint knut_dgeevx (char balanc, char jobvl, char jobvr, char sense, ptrdiff_t n,
	     double *a, ptrdiff_t lda, double *wr, double *wi, double *vl, ptrdiff_t ldvl,
	     double *vr, ptrdiff_t ldvr, blasint *ilo, blasint *ihi, double *scale,
	     double *abnrm, double *rconde, double *rcondv, double *work,
	     blasint *lwork, blasint *iwork, blasint *info)
{
  blasint n_ = (blasint)n;
  blasint lda_ = (blasint)lda;
  blasint ldvl_ = (blasint)ldvl;
  blasint ldvr_ = (blasint)ldvr;
  return dgeevx_(&balanc, &jobvl, &jobvr, &sense, &n_,
    a, &lda_, wr, wi, vl, &ldvl_,
    vr, &ldvr_, ilo, ihi, scale,
    abnrm, rconde, rcondv, work,
    lwork, iwork, info, 1, 1, 1, 1);
}

inline blasint knut_dgeev (char jobvl, char jobvr, ptrdiff_t n, double *a, ptrdiff_t lda, double *wr,
	    double *wi, double *vl, ptrdiff_t ldvl, double *vr, ptrdiff_t ldvr,
	    double *work, ptrdiff_t lwork, blasint *info)
{
  blasint n_ = (blasint)n;
  blasint lda_ = (blasint)lda;
  blasint ldvl_ = (blasint)ldvl;
  blasint ldvr_ = (blasint)ldvr;
  blasint lwork_ = (blasint)lwork;
  return dgeev_(&jobvl, &jobvr, &n_, a, &lda_, wr, 
    wi, vl, &ldvl_, vr, &ldvr_,
    work, &lwork_, info, 1, 1);
}

inline blasint knut_dlarnv (blasint *idist, blasint *iseed, blasint *n, double *x)
{
  return dlarnv_(idist, iseed, n, x);
}

inline double knut_dlamch (const char *cmach, blasint cmach_len)
{
  return dlamch_(cmach, cmach_len);
}

#endif
