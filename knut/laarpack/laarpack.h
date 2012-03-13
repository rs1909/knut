#ifndef LAARPACK_H
#define LAARPACK_H

int dneupd_(int * rvec, char *howmny, bool * select, double *dr, double *di,
	     double *z__, int *ldz, double *sigmar, double *sigmai,
	     double *workev, char *bmat, int *n, char *which, int *nev,
	     double *tol, double *resid, int *ncv, double *v, int *ldv,
	     int *iparam, int *ipntr, double *workd, double *workl, int *lworkl,
	     int *info, int howmny_len, int bmat_len, int which_len);

int dnaupd_(int *ido, char *bmat, int *n, char *which, int *nev, double *tol,
	     double *resid, int *ncv, double *v, int *ldv, int *iparam,
	     int *ipntr, double *workd, double *workl, int *lworkl, int *info,
	     int bmat_len, int which_len);

int dgesvx_(char *fact, char *trans, int *n, int *nrhs, double *a, int *lda,
	     double *af, int *ldaf, int *ipiv, char *equed, double *r__,
	     double *c__, double *b, int *ldb, double *x, int *ldx,
	     double *rcond, double *ferr, double *berr, double *work,
	     int *iwork, int *info, int fact_len, int trans_len, int equed_len);

int dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b,
	    int *ldb, int *info);

int dgeevx_(char *balanc, char *jobvl, char *jobvr, char *sense, int *n,
	     double *a, int *lda, double *wr, double *wi, double *vl, int *ldvl,
	     double *vr, int *ldvr, int *ilo, int *ihi, double *scale,
	     double *abnrm, double *rconde, double *rcondv, double *work,
	     int *lwork, int *iwork, int *info, int balanc_len, int jobvl_len,
	     int jobvr_len, int sense_len);

int dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda, double *wr,
	    double *wi, double *vl, int *ldvl, double *vr, int *ldvr,
	    double *work, int *lwork, int *info, int jobvl_len, int jobvr_len);

int dlarnv_(int *idist, int *iseed, int *n, double *x);

double dlamch_(const char *cmach, int cmach_len);

// for matrix.cpp and spmatrix.cpp

inline int knut_dneupd (bool rvec, char howmny, bool * select, double *dr, double *di,
	     double *z__, size_t ldz, double *sigmar, double *sigmai,
	     double *workev, char bmat, size_t n, char which1, char which2, size_t nev,
	     double *tol, double *resid, size_t ncv, double *v, size_t ldv,
	     int *iparam, int *ipntr, double *workd, double *workl, size_t lworkl,
	     int *info)
{
  int n_ = (int)n;
  int rvec_ = (int)rvec;
  int ldz_ = (int)ldz;
  int nev_ = (int)nev;  
  int ncv_ = (int)ncv;  
  int ldv_ = (int)ldv;
  int lworkl_ = (int)lworkl;
  char which[2] = {which1, which2};  
  return dneupd_(&rvec_, &howmny, select, dr, di,
    z__, &ldz_, sigmar, sigmai,
    workev, &bmat, &n_, which, &nev_,
    tol, resid, &ncv_, v, &ldv_,
    iparam, ipntr, workd, workl, &lworkl_,
    info, 1, 1, 2);
}

inline int knut_dnaupd (int *ido, char bmat, size_t n, char which1, char which2, size_t nev, double *tol,
	     double *resid, size_t ncv, double *v, size_t ldv, int *iparam,
	     int *ipntr, double *workd, double *workl, size_t lworkl, int *info)
{
  int n_ = (int)n;
  int nev_ = (int)nev;  
  int ncv_ = (int)ncv;  
  int ldv_ = (int)ldv;
  int lworkl_ = (int)lworkl;
  char which[2] = {which1, which2}; 
  return dnaupd_(ido, &bmat, &n_, which, &nev_, tol,
    resid, &ncv_, v, &ldv_, iparam,
    ipntr, workd, workl, &lworkl_, info,
    1, 2);
}

inline int knut_dgesvx (char fact, char trans, size_t n, size_t nrhs, double *a, size_t lda,
	     double *af, size_t ldaf, int *ipiv, char equed, double *r__,
	     double *c__, double *b, size_t ldb, double *x, size_t ldx,
	     double *rcond, double *ferr, double *berr, double *work,
	     int *iwork, int *info)
{
  int n_ = (int)n;
  int nrhs_ = (int)nrhs;
  int lda_ = (int)lda;
  int ldb_ = (int)ldb;
  int ldx_ = (int)ldx;
  int ldaf_ = (int)ldaf;
  return dgesvx_(&fact, &trans, &n_, &nrhs_, a, &lda_,
    af, &ldaf_, ipiv, &equed, r__,
    c__, b, &ldb_, x, &ldx_,
    rcond, ferr, berr, work,
    iwork, info, 1, 1, 1);
}

inline int knut_dgesv (int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b,
	    int *ldb, int *info)
{
  return dgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
}

inline int knut_dgeevx (char balanc, char jobvl, char jobvr, char sense, size_t n,
	     double *a, size_t lda, double *wr, double *wi, double *vl, size_t ldvl,
	     double *vr, size_t ldvr, int *ilo, int *ihi, double *scale,
	     double *abnrm, double *rconde, double *rcondv, double *work,
	     int *lwork, int *iwork, int *info)
{
  int n_ = (int)n;
  int lda_ = (int)lda;
  int ldvl_ = (int)ldvl;
  int ldvr_ = (int)ldvr;
  return dgeevx_(&balanc, &jobvl, &jobvr, &sense, &n_,
    a, &lda_, wr, wi, vl, &ldvl_,
    vr, &ldvr_, ilo, ihi, scale,
    abnrm, rconde, rcondv, work,
    lwork, iwork, info, 1, 1, 1, 1);
}

inline int knut_dgeev (char jobvl, char jobvr, size_t n, double *a, size_t lda, double *wr,
	    double *wi, double *vl, size_t ldvl, double *vr, size_t ldvr,
	    double *work, size_t lwork, int *info)
{
  int n_ = (int)n;
  int lda_ = (int)lda;
  int ldvl_ = (int)ldvl;
  int ldvr_ = (int)ldvr;
  int lwork_ = (int)lwork;
  return dgeev_(&jobvl, &jobvr, &n_, a, &lda_, wr, 
    wi, vl, &ldvl_, vr, &ldvr_,
    work, &lwork_, info, 1, 1);
}

inline int knut_dlarnv (int *idist, int *iseed, int *n, double *x)
{
  return dlarnv_(idist, iseed, n, x);
}

inline double knut_dlamch (const char *cmach, int cmach_len)
{
  return dlamch_(cmach, cmach_len);
}

#endif
