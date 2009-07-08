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

inline int knut_dneupd (int * rvec, char *howmny, bool * select, double *dr, double *di,
	     double *z__, int *ldz, double *sigmar, double *sigmai,
	     double *workev, char *bmat, int *n, char *which, int *nev,
	     double *tol, double *resid, int *ncv, double *v, int *ldv,
	     int *iparam, int *ipntr, double *workd, double *workl, int *lworkl,
	     int *info, int howmny_len, int bmat_len, int which_len)
{
  return dneupd_(rvec, howmny, select, dr, di,
    z__, ldz, sigmar, sigmai,
    workev, bmat, n, which, nev,
    tol, resid, ncv, v, ldv,
    iparam, ipntr, workd, workl, lworkl,
    info, howmny_len, bmat_len, which_len);
}

inline int knut_dnaupd (int *ido, char *bmat, int *n, char *which, int *nev, double *tol,
	     double *resid, int *ncv, double *v, int *ldv, int *iparam,
	     int *ipntr, double *workd, double *workl, int *lworkl, int *info,
	     int bmat_len, int which_len)
{
  return dnaupd_(ido, bmat, n, which, nev, tol,
    resid, ncv, v, ldv, iparam,
    ipntr, workd, workl, lworkl, info,
    bmat_len, which_len);
}

inline int knut_dgesvx (char *fact, char *trans, int *n, int *nrhs, double *a, int *lda,
	     double *af, int *ldaf, int *ipiv, char *equed, double *r__,
	     double *c__, double *b, int *ldb, double *x, int *ldx,
	     double *rcond, double *ferr, double *berr, double *work,
	     int *iwork, int *info, int fact_len, int trans_len, int equed_len)
{
  return dgesvx_(fact, trans, n, nrhs, a, lda,
    af, ldaf, ipiv, equed, r__,
    c__, b, ldb, x, ldx,
    rcond, ferr, berr, work,
    iwork, info, fact_len, trans_len, equed_len);
}

inline int knut_dgesv (int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b,
	    int *ldb, int *info)
{
  return dgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
}

inline int knut_dgeevx (char *balanc, char *jobvl, char *jobvr, char *sense, int *n,
	     double *a, int *lda, double *wr, double *wi, double *vl, int *ldvl,
	     double *vr, int *ldvr, int *ilo, int *ihi, double *scale,
	     double *abnrm, double *rconde, double *rcondv, double *work,
	     int *lwork, int *iwork, int *info, int balanc_len, int jobvl_len,
	     int jobvr_len, int sense_len)
{
  return dgeevx_(balanc, jobvl, jobvr, sense, n,
    a, lda, wr, wi, vl, ldvl,
    vr, ldvr, ilo, ihi, scale,
    abnrm, rconde, rcondv, work,
    lwork, iwork, info, balanc_len, jobvl_len,
    jobvr_len, sense_len);
}

inline int knut_dgeev (char *jobvl, char *jobvr, int *n, double *a, int *lda, double *wr,
	    double *wi, double *vl, int *ldvl, double *vr, int *ldvr,
	    double *work, int *lwork, int *info, int jobvl_len, int jobvr_len)
{
  return dgeev_(jobvl, jobvr, n, a, lda, wr, 
    wi, vl, ldvl, vr, ldvr,
    work, lwork, info, jobvl_len, jobvr_len);
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
