#ifndef LAARPACK_H
#define LAARPACK_H

int knut_dneupd (bool * rvec, char *howmny, bool * select, double *dr, double *di,
	     double *z__, int *ldz, double *sigmar, double *sigmai,
	     double *workev, char *bmat, int *n, char *which, int *nev,
	     double *tol, double *resid, int *ncv, double *v, int *ldv,
	     int *iparam, int *ipntr, double *workd, double *workl, int *lworkl,
	     int *info, int howmny_len, int bmat_len, int which_len);

int knut_dnaupd (int *ido, char *bmat, int *n, char *which, int *nev, double *tol,
	     double *resid, int *ncv, double *v, int *ldv, int *iparam,
	     int *ipntr, double *workd, double *workl, int *lworkl, int *info,
	     int bmat_len, int which_len);

int knut_dgesvx (char *fact, char *trans, int *n, int *nrhs, double *a, int *lda,
	     double *af, int *ldaf, int *ipiv, char *equed, double *r__,
	     double *c__, double *b, int *ldb, double *x, int *ldx,
	     double *rcond, double *ferr, double *berr, double *work,
	     int *iwork, int *info, int fact_len, int trans_len, int equed_len);

int knut_dgesv (int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b,
	    int *ldb, int *info);

int knut_dgeevx (char *balanc, char *jobvl, char *jobvr, char *sense, int *n,
	     double *a, int *lda, double *wr, double *wi, double *vl, int *ldvl,
	     double *vr, int *ldvr, int *ilo, int *ihi, double *scale,
	     double *abnrm, double *rconde, double *rcondv, double *work,
	     int *lwork, int *iwork, int *info, int balanc_len, int jobvl_len,
	     int jobvr_len, int sense_len);

int knut_dgeev (char *jobvl, char *jobvr, int *n, double *a, int *lda, double *wr,
	    double *wi, double *vl, int *ldvl, double *vr, int *ldvr,
	    double *work, int *lwork, int *info, int jobvl_len, int jobvr_len);

int knut_dlarnv (int *idist, int *iseed, int *n, double *x);

double knut_dlamch (const char *cmach, int cmach_len);

#endif
