/* ARPACK/laarpack.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
   on Microsoft Windows system, link with libf2c.lib;
   on Linux or Unix systems, link with .../path/to/libf2c.a -lm
   or, if you install libf2c.a in a standard place, with -lf2c -lm
   -- in that order, at the end of the command line, as in
   cc *.o -lf2c -lm
   Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

   http://www.netlib.org/f2c/libf2c.zip
 */

#include "f2c.h"
#include "cblaswrap.h"

int dneupd_ (logical * rvec, char *howmny, logical * select, doublereal * dr, doublereal * di, doublereal * z__, integer * ldz, doublereal * sigmar,
				 doublereal * sigmai, doublereal * workev, char *bmat, integer * n, char *which, integer * nev, doublereal * tol, doublereal * resid, integer * ncv,
				 doublereal * v, integer * ldv, integer * iparam, integer * ipntr, doublereal * workd, doublereal * workl, integer * lworkl, integer * info,
				 ftnlen howmny_len, ftnlen bmat_len, ftnlen which_len);

int dnaupd_ (integer * ido, char *bmat, integer * n, char *which, integer * nev, doublereal * tol, doublereal * resid, integer * ncv, doublereal * v,
				 integer * ldv, integer * iparam, integer * ipntr, doublereal * workd, doublereal * workl, integer * lworkl, integer * info, ftnlen bmat_len,
				 ftnlen which_len);

int dgesvx_ (char *fact, char *trans, integer * n, integer * nrhs, doublereal * a, integer * lda, doublereal * af, integer * ldaf, integer * ipiv, char *equed,
				 doublereal * r__, doublereal * c__, doublereal * b, integer * ldb, doublereal * x, integer * ldx, doublereal * rcond, doublereal * ferr,
				 doublereal * berr, doublereal * work, integer * iwork, integer * info, ftnlen fact_len, ftnlen trans_len, ftnlen equed_len);

int dgesv_ (integer * n, integer * nrhs, doublereal * a, integer * lda, integer * ipiv, doublereal * b, integer * ldb, integer * info);

int dgeevx_ (char *balanc, char *jobvl, char *jobvr, char *sense, integer * n, doublereal * a, integer * lda, doublereal * wr, doublereal * wi, doublereal * vl,
				 integer * ldvl, doublereal * vr, integer * ldvr, integer * ilo, integer * ihi, doublereal * scale, doublereal * abnrm, doublereal * rconde,
				 doublereal * rcondv, doublereal * work, integer * lwork, integer * iwork, integer * info, ftnlen balanc_len, ftnlen jobvl_len, ftnlen jobvr_len,
				 ftnlen sense_len);

int dgeev_ (char *jobvl, char *jobvr, integer * n, doublereal * a, integer * lda, doublereal * wr, doublereal * wi, doublereal * vl, integer * ldvl,
				doublereal * vr, integer * ldvr, doublereal * work, integer * lwork, integer * info, ftnlen jobvl_len, ftnlen jobvr_len);

// from dlamch.c
doublereal dlamch_ (char *cmach, ftnlen cmach_len);

// static
static int dstatn_ (void);

static int dsortc_ (char *which, logical * apply, integer * n, doublereal * xreal, doublereal * ximag, doublereal * y, ftnlen which_len);

static int dngets_ (integer * ishift, char *which, integer * kev, integer * np, doublereal * ritzr, doublereal * ritzi, doublereal * bounds,
						  doublereal * shiftr, doublereal * shifti, ftnlen which_len);

static int dnconv_ (integer * n, doublereal * ritzr, doublereal * ritzi, doublereal * bounds, doublereal * tol, integer * nconv);

static int dneigh_ (doublereal * rnorm, integer * n, doublereal * h__, integer * ldh, doublereal * ritzr, doublereal * ritzi, doublereal * bounds,
						  doublereal * q, integer * ldq, doublereal * workl, integer * ierr);

static int dnaup2_ (integer * ido, char *bmat, integer * n, char *which, integer * nev, integer * np, doublereal * tol, doublereal * resid, integer * mode,
						  integer * iupd, integer * ishift, integer * mxiter, doublereal * v, integer * ldv, doublereal * h__, integer * ldh, doublereal * ritzr,
						  doublereal * ritzi, doublereal * bounds, doublereal * q, integer * ldq, doublereal * workl, integer * ipntr, doublereal * workd,
						  integer * info, ftnlen bmat_len, ftnlen which_len);

static int dnapps_ (integer * n, integer * kev, integer * np, doublereal * shiftr, doublereal * shifti, doublereal * v, integer * ldv, doublereal * h__,
						  integer * ldh, doublereal * resid, doublereal * q, integer * ldq, doublereal * workl, doublereal * workd);

static int ivout_ (integer * lout, integer * n, integer * ix, integer * idigit, char *ifmt, ftnlen ifmt_len);

static int dnaitr_ (integer * ido, char *bmat, integer * n, integer * k, integer * np, integer * nb, doublereal * resid, doublereal * rnorm, doublereal * v,
						  integer * ldv, doublereal * h__, integer * ldh, integer * ipntr, doublereal * workd, integer * info, ftnlen bmat_len);

static int dmout_ (integer * lout, integer * m, integer * n, doublereal * a, integer * lda, integer * idigit, char *ifmt, ftnlen ifmt_len);

static int dlaqrb_ (logical * wantt, integer * n, integer * ilo, integer * ihi, doublereal * h__, integer * ldh, doublereal * wr, doublereal * wi,
						  doublereal * z__, integer * info);

static int second_ (real * t);

static int dvout_ (integer * lout, integer * n, doublereal * sx, integer * idigit, char *ifmt, ftnlen ifmt_len);

static int dgetv0_ (integer * ido, char *bmat, integer * itry, logical * initv, integer * n, integer * j, doublereal * v, integer * ldv, doublereal * resid,
						  doublereal * rnorm, integer * ipntr, doublereal * workd, integer * ierr, ftnlen bmat_len);

static integer ieeeck_ (integer * ispec, real * zero, real * one);

static int dtrsyl_ (char *trana, char *tranb, integer * isgn, integer * m, integer * n, doublereal * a, integer * lda, doublereal * b, integer * ldb,
						  doublereal * c__, integer * ldc, doublereal * scale, integer * info, ftnlen trana_len, ftnlen tranb_len);

static int dtrsen_ (char *job, char *compq, logical * select, integer * n, doublereal * t, integer * ldt, doublereal * q, integer * ldq, doublereal * wr,
						  doublereal * wi, integer * m, doublereal * s, doublereal * sep, doublereal * work, integer * lwork, integer * iwork, integer * liwork,
						  integer * info, ftnlen job_len, ftnlen compq_len);

static int dtrexc_ (char *compq, integer * n, doublereal * t, integer * ldt, doublereal * q, integer * ldq, integer * ifst, integer * ilst, doublereal * work,
						  integer * info, ftnlen compq_len);

static int dorm2r_ (char *side, char *trans, integer * m, integer * n, integer * k, doublereal * a, integer * lda, doublereal * tau, doublereal * c__,
						  integer * ldc, doublereal * work, integer * info, ftnlen side_len, ftnlen trans_len);

static int dorgqr_ (integer * m, integer * n, integer * k, doublereal * a, integer * lda, doublereal * tau, doublereal * work, integer * lwork, integer * info);

static int dorg2r_ (integer * m, integer * n, integer * k, doublereal * a, integer * lda, doublereal * tau, doublereal * work, integer * info);

static int dlaruv_ (integer * iseed, integer * n, doublereal * x);

static int dlarnv_ (integer * idist, integer * iseed, integer * n, doublereal * x);

static int dlarft_ (char *direct, char *storev, integer * n, integer * k, doublereal * v, integer * ldv, doublereal * tau, doublereal * t, integer * ldt,
						  ftnlen direct_len, ftnlen storev_len);

static int dlaqtr_ (logical * ltran, logical * lreal, integer * n, doublereal * t, integer * ldt, doublereal * b, doublereal * w, doublereal * scale,
						  doublereal * x, doublereal * work, integer * info);

static int dlassq_ (integer * n, doublereal * x, integer * incx, doublereal * scale, doublereal * sumsq);

static int dlaln2_ (logical * ltrans, integer * na, integer * nw, doublereal * smin, doublereal * ca, doublereal * a, integer * lda, doublereal * d1,
						  doublereal * d2, doublereal * b, integer * ldb, doublereal * wr, doublereal * wi, doublereal * x, integer * ldx, doublereal * scale,
						  doublereal * xnorm, integer * info);

static int dlasy2_ (logical * ltranl, logical * ltranr, integer * isgn, integer * n1, integer * n2, doublereal * tl, integer * ldtl, doublereal * tr,
						  integer * ldtr, doublereal * b, integer * ldb, doublereal * scale, doublereal * x, integer * ldx, doublereal * xnorm, integer * info);

static int dlanv2_ (doublereal * a, doublereal * b, doublereal * c__, doublereal * d__, doublereal * rt1r, doublereal * rt1i, doublereal * rt2r,
						  doublereal * rt2i, doublereal * cs, doublereal * sn);

static int dlaexc_ (logical * wantq, integer * n, doublereal * t, integer * ldt, doublereal * q, integer * ldq, integer * j1, integer * n1, integer * n2,
						  doublereal * work, integer * info);

static int dladiv_ (doublereal * a, doublereal * b, doublereal * c__, doublereal * d__, doublereal * p, doublereal * q);

static int dlarfx_ (char *side, integer * m, integer * n, doublereal * v, doublereal * tau, doublereal * c__, integer * ldc, doublereal * work,
						  ftnlen side_len);

static int dlaset_ (char *uplo, integer * m, integer * n, doublereal * alpha, doublereal * beta, doublereal * a, integer * lda, ftnlen uplo_len);

static int dlahqr_ (logical * wantt, logical * wantz, integer * n, integer * ilo, integer * ihi, doublereal * h__, integer * ldh, doublereal * wr,
						  doublereal * wi, integer * iloz, integer * ihiz, doublereal * z__, integer * ldz, integer * info);

static doublereal dlanhs_ (char *norm, integer * n, doublereal * a, integer * lda, doublereal * work, ftnlen norm_len);

static int dlaswp_ (integer * n, doublereal * a, integer * lda, integer * k1, integer * k2, integer * ipiv, integer * incx);

static int dgetf2_ (integer * m, integer * n, doublereal * a, integer * lda, integer * ipiv, integer * info);

static doublereal dlantr_ (char *norm, char *uplo, char *diag, integer * m, integer * n, doublereal * a, integer * lda, doublereal * work, ftnlen norm_len,
									ftnlen uplo_len, ftnlen diag_len);

static int dlaqge_ (integer * m, integer * n, doublereal * a, integer * lda, doublereal * r__, doublereal * c__, doublereal * rowcnd, doublereal * colcnd,
						  doublereal * amax, char *equed, ftnlen equed_len);

static int dgetrf_ (integer * m, integer * n, doublereal * a, integer * lda, integer * ipiv, integer * info);

static int dgetrs_ (char *trans, integer * n, integer * nrhs, doublereal * a, integer * lda, integer * ipiv, doublereal * b, integer * ldb, integer * info,
						  ftnlen trans_len);

static int dgerfs_ (char *trans, integer * n, integer * nrhs, doublereal * a, integer * lda, doublereal * af, integer * ldaf, integer * ipiv, doublereal * b,
						  integer * ldb, doublereal * x, integer * ldx, doublereal * ferr, doublereal * berr, doublereal * work, integer * iwork, integer * info,
						  ftnlen trans_len);

static int dgeqr2_ (integer * m, integer * n, doublereal * a, integer * lda, doublereal * tau, doublereal * work, integer * info);

static int dlahrd_ (integer * n, integer * k, integer * nb, doublereal * a, integer * lda, doublereal * tau, doublereal * t, integer * ldt, doublereal * y,
						  integer * ldy);

static int dlarfb_ (char *side, char *trans, char *direct, char *storev, integer * m, integer * n, integer * k, doublereal * v, integer * ldv, doublereal * t,
						  integer * ldt, doublereal * c__, integer * ldc, doublereal * work, integer * ldwork, ftnlen side_len, ftnlen trans_len, ftnlen direct_len,
						  ftnlen storev_len);

static int dlarfg_ (integer * n, doublereal * alpha, doublereal * x, integer * incx, doublereal * tau);

static int dlarf_ (char *side, integer * m, integer * n, doublereal * v, integer * incv, doublereal * tau, doublereal * c__, integer * ldc, doublereal * work,
						 ftnlen side_len);

static int dgehd2_ (integer * n, integer * ilo, integer * ihi, doublereal * a, integer * lda, doublereal * tau, doublereal * work, integer * info);

static int dtrsna_ (char *job, char *howmny, logical * select, integer * n, doublereal * t, integer * ldt, doublereal * vl, integer * ldvl, doublereal * vr,
						  integer * ldvr, doublereal * s, doublereal * sep, integer * mm, integer * m, doublereal * work, integer * ldwork, integer * iwork,
						  integer * info, ftnlen job_len, ftnlen howmny_len);

static int dtrevc_ (char *side, char *howmny, logical * select, integer * n, doublereal * t, integer * ldt, doublereal * vl, integer * ldvl, doublereal * vr,
						  integer * ldvr, integer * mm, integer * m, doublereal * work, integer * info, ftnlen side_len, ftnlen howmny_len);

static int dhseqr_ (char *job, char *compz, integer * n, integer * ilo, integer * ihi, doublereal * h__, integer * ldh, doublereal * wr, doublereal * wi,
						  doublereal * z__, integer * ldz, doublereal * work, integer * lwork, integer * info, ftnlen job_len, ftnlen compz_len);

static int dorghr_ (integer * n, integer * ilo, integer * ihi, doublereal * a, integer * lda, doublereal * tau, doublereal * work, integer * lwork,
						  integer * info);

static integer ilaenv_ (integer * ispec, char *name__, char *opts, integer * n1, integer * n2, integer * n3, integer * n4, ftnlen name_len, ftnlen opts_len);

static int dlartg_ (doublereal * f, doublereal * g, doublereal * cs, doublereal * sn, doublereal * r__);

static int dlacpy_ (char *uplo, integer * m, integer * n, doublereal * a, integer * lda, doublereal * b, integer * ldb, ftnlen uplo_len);

static int dlascl_ (char *type__, integer * kl, integer * ku, doublereal * cfrom, doublereal * cto, integer * m, integer * n, doublereal * a, integer * lda,
						  integer * info, ftnlen type_len);

static int dgehrd_ (integer * n, integer * ilo, integer * ihi, doublereal * a, integer * lda, doublereal * tau, doublereal * work, integer * lwork,
						  integer * info);

static doublereal dlange_ (char *norm, integer * m, integer * n, doublereal * a, integer * lda, doublereal * work, ftnlen norm_len);

static int dlabad_ (doublereal * small, doublereal * large);

static doublereal dlapy2_ (doublereal * x, doublereal * y);

static int dgeequ_ (integer * m, integer * n, doublereal * a, integer * lda, doublereal * r__, doublereal * c__, doublereal * rowcnd, doublereal * colcnd,
						  doublereal * amax, integer * info);

static int dlatrs_ (char *uplo, char *trans, char *diag, char *normin, integer * n, doublereal * a, integer * lda, doublereal * x, doublereal * scale,
						  doublereal * cnorm, integer * info, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len, ftnlen normin_len);

static int dlacon_ (integer * n, doublereal * v, doublereal * x, integer * isgn, doublereal * est, integer * kase);

static int drscl_ (integer * n, doublereal * sa, doublereal * sx, integer * incx);

static int dgecon_ (char *norm, integer * n, doublereal * a, integer * lda, doublereal * anorm, doublereal * rcond, doublereal * work, integer * iwork,
						  integer * info, ftnlen norm_len);

static int dgebal_ (char *job, integer * n, doublereal * a, integer * lda, integer * ilo, integer * ihi, doublereal * scale, integer * info, ftnlen job_len);

static int xerbla_ (char *srname, integer * info, ftnlen srname_len);

static int dgebak_ (char *job, char *side, integer * n, integer * ilo, integer * ihi, doublereal * scale, integer * m, doublereal * v, integer * ldv,
						  integer * info, ftnlen job_len, ftnlen side_len);

doublereal etime_ (real *);

/* Common Block Declarations */

struct
{
	integer logfil, ndigit, mgetv0, msaupd, msaup2, msaitr, mseigt, msapps,
		msgets, mseupd, mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd, mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd;
} debug_;

#define debug_1 debug_

struct
{
	integer nopx, nbx, nrorth, nitref, nrstrt;
	real tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv, tnaupd,
		tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv, tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv, tmvopx, tmvbx, tgetv0, titref, trvec;
} timing_;

#define timing_1 timing_

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c__8 = 8;
static integer c_n1 = -1;
static integer c__4 = 4;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c__65 = 65;
static doublereal c_b347 = -1.;
static doublereal c_b348 = 1.;
static doublereal c_b507 = 0.;
static integer c__15 = 15;
static logical c_false = FALSE_;
static logical c_true = TRUE_;
static integer c__16 = 16;
static doublereal c_b1496 = .5;
static real c_b2255 = 0.f;
static real c_b2256 = 1.f;
static doublereal c_b2616 = .66666666666666663;

static int
dgebak_ (char *job, char *side, integer * n, integer * ilo,
			integer * ihi, doublereal * scale, integer * m, doublereal * v, integer * ldv, integer * info, ftnlen job_len, ftnlen side_len)
{
	/* System generated locals */
	integer v_dim1, v_offset, i__1;

	/* Local variables */
	static integer i__, k;
	static doublereal s;
	static integer ii;

	static logical leftv;

	static logical rightv;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     September 30, 1994 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DGEBAK forms the right or left eigenvectors of a real general matrix */
	/*  by backward transformation on the computed eigenvectors of the */
	/*  balanced matrix output by DGEBAL. */

	/*  Arguments */
	/*  ========= */

	/*  JOB     (input) CHARACTER*1 */
	/*          Specifies the type of backward transformation required: */
	/*          = 'N', do nothing, return immediately; */
	/*          = 'P', do backward transformation for permutation only; */
	/*          = 'S', do backward transformation for scaling only; */
	/*          = 'B', do backward transformations for both permutation and */
	/*                 scaling. */
	/*          JOB must be the same as the argument JOB supplied to DGEBAL. */

	/*  SIDE    (input) CHARACTER*1 */
	/*          = 'R':  V contains right eigenvectors; */
	/*          = 'L':  V contains left eigenvectors. */

	/*  N       (input) INTEGER */
	/*          The number of rows of the matrix V.  N >= 0. */

	/*  ILO     (input) INTEGER */
	/*  IHI     (input) INTEGER */
	/*          The integers ILO and IHI determined by DGEBAL. */
	/*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */

	/*  SCALE   (input) DOUBLE PRECISION array, dimension (N) */
	/*          Details of the permutation and scaling factors, as returned */
	/*          by DGEBAL. */

	/*  M       (input) INTEGER */
	/*          The number of columns of the matrix V.  M >= 0. */

	/*  V       (input/output) DOUBLE PRECISION array, dimension (LDV,M) */
	/*          On entry, the matrix of right or left eigenvectors to be */
	/*          transformed, as returned by DHSEIN or DTREVC. */
	/*          On exit, V is overwritten by the transformed eigenvectors. */

	/*  LDV     (input) INTEGER */
	/*          The leading dimension of the array V. LDV >= max(1,N). */

	/*  INFO    (output) INTEGER */
	/*          = 0:  successful exit */
	/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Decode and Test the input parameters */

	/* Parameter adjustments */
	--scale;
	v_dim1 = *ldv;
	v_offset = 1 + v_dim1;
	v -= v_offset;

	/* Function Body */
	rightv = lsame_ (side, "R", (ftnlen) 1, (ftnlen) 1);
	leftv = lsame_ (side, "L", (ftnlen) 1, (ftnlen) 1);

	*info = 0;
	if (!lsame_ (job, "N", (ftnlen) 1, (ftnlen) 1) && !lsame_ (job, "P", (ftnlen) 1, (ftnlen) 1) && !lsame_ (job, "S", (ftnlen) 1, (ftnlen) 1)
		 && !lsame_ (job, "B", (ftnlen) 1, (ftnlen) 1))
	{
		*info = -1;
	}
	else if (!rightv && !leftv)
	{
		*info = -2;
	}
	else if (*n < 0)
	{
		*info = -3;
	}
	else if (*ilo < 1 || *ilo > max (1, *n))
	{
		*info = -4;
	}
	else if (*ihi < min (*ilo, *n) || *ihi > *n)
	{
		*info = -5;
	}
	else if (*m < 0)
	{
		*info = -7;
	}
	else if (*ldv < max (1, *n))
	{
		*info = -9;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DGEBAK", &i__1, (ftnlen) 6);
		return 0;
	}

	/*     Quick return if possible */

	if (*n == 0)
	{
		return 0;
	}
	if (*m == 0)
	{
		return 0;
	}
	if (lsame_ (job, "N", (ftnlen) 1, (ftnlen) 1))
	{
		return 0;
	}

	if (*ilo == *ihi)
	{
		goto L30;
	}

	/*     Backward balance */

	if (lsame_ (job, "S", (ftnlen) 1, (ftnlen) 1) || lsame_ (job, "B", (ftnlen) 1, (ftnlen) 1))
	{

		if (rightv)
		{
			i__1 = *ihi;
			for (i__ = *ilo; i__ <= i__1; ++i__)
			{
				s = scale[i__];
				dscal_ (m, &s, &v[i__ + v_dim1], ldv);
				/* L10: */
			}
		}

		if (leftv)
		{
			i__1 = *ihi;
			for (i__ = *ilo; i__ <= i__1; ++i__)
			{
				s = 1. / scale[i__];
				dscal_ (m, &s, &v[i__ + v_dim1], ldv);
				/* L20: */
			}
		}

	}

	/*     Backward permutation */

	/*     For  I = ILO-1 step -1 until 1, */
	/*              IHI+1 step 1 until N do -- */

 L30:
	if (lsame_ (job, "P", (ftnlen) 1, (ftnlen) 1) || lsame_ (job, "B", (ftnlen) 1, (ftnlen) 1))
	{
		if (rightv)
		{
			i__1 = *n;
			for (ii = 1; ii <= i__1; ++ii)
			{
				i__ = ii;
				if (i__ >= *ilo && i__ <= *ihi)
				{
					goto L40;
				}
				if (i__ < *ilo)
				{
					i__ = *ilo - ii;
				}
				k = (integer) scale[i__];
				if (k == i__)
				{
					goto L40;
				}
				dswap_ (m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
			 L40:
				;
			}
		}

		if (leftv)
		{
			i__1 = *n;
			for (ii = 1; ii <= i__1; ++ii)
			{
				i__ = ii;
				if (i__ >= *ilo && i__ <= *ihi)
				{
					goto L50;
				}
				if (i__ < *ilo)
				{
					i__ = *ilo - ii;
				}
				k = (integer) scale[i__];
				if (k == i__)
				{
					goto L50;
				}
				dswap_ (m, &v[i__ + v_dim1], ldv, &v[k + v_dim1], ldv);
			 L50:
				;
			}
		}
	}

	return 0;

	/*     End of DGEBAK */

}	/* dgebak_ */

static int
dgebal_ (char *job, integer * n, doublereal * a, integer * lda, integer * ilo, integer * ihi, doublereal * scale, integer * info, ftnlen job_len)
{
	/* System generated locals */
	integer a_dim1, a_offset, i__1, i__2;
	doublereal d__1, d__2;

	/* Local variables */
	static doublereal c__, f, g;
	static integer i__, j, k, l, m;
	static doublereal r__, s, ca, ra;
	static integer ica, ira, iexc;

	static doublereal sfmin1, sfmin2, sfmax1, sfmax2;
	static logical noconv;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     June 30, 1999 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DGEBAL balances a general real matrix A.  This involves, first, */
	/*  permuting A by a similarity transformation to isolate eigenvalues */
	/*  in the first 1 to ILO-1 and last IHI+1 to N elements on the */
	/*  diagonal; and second, applying a diagonal similarity transformation */
	/*  to rows and columns ILO to IHI to make the rows and columns as */
	/*  close in norm as possible.  Both steps are optional. */

	/*  Balancing may reduce the 1-norm of the matrix, and improve the */
	/*  accuracy of the computed eigenvalues and/or eigenvectors. */

	/*  Arguments */
	/*  ========= */

	/*  JOB     (input) CHARACTER*1 */
	/*          Specifies the operations to be performed on A: */
	/*          = 'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0 */
	/*                  for i = 1,...,N; */
	/*          = 'P':  permute only; */
	/*          = 'S':  scale only; */
	/*          = 'B':  both permute and scale. */

	/*  N       (input) INTEGER */
	/*          The order of the matrix A.  N >= 0. */

	/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          On entry, the input matrix A. */
	/*          On exit,  A is overwritten by the balanced matrix. */
	/*          If JOB = 'N', A is not referenced. */
	/*          See Further Details. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(1,N). */

	/*  ILO     (output) INTEGER */
	/*  IHI     (output) INTEGER */
	/*          ILO and IHI are set to integers such that on exit */
	/*          A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N. */
	/*          If JOB = 'N' or 'S', ILO = 1 and IHI = N. */

	/*  SCALE   (output) DOUBLE PRECISION array, dimension (N) */
	/*          Details of the permutations and scaling factors applied to */
	/*          A.  If P(j) is the index of the row and column interchanged */
	/*          with row and column j and D(j) is the scaling factor */
	/*          applied to row and column j, then */
	/*          SCALE(j) = P(j)    for j = 1,...,ILO-1 */
	/*                   = D(j)    for j = ILO,...,IHI */
	/*                   = P(j)    for j = IHI+1,...,N. */
	/*          The order in which the interchanges are made is N to IHI+1, */
	/*          then 1 to ILO-1. */

	/*  INFO    (output) INTEGER */
	/*          = 0:  successful exit. */
	/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */

	/*  Further Details */
	/*  =============== */

	/*  The permutations consist of row and column interchanges which put */
	/*  the matrix in the form */

	/*             ( T1   X   Y  ) */
	/*     P A P = (  0   B   Z  ) */
	/*             (  0   0   T2 ) */

	/*  where T1 and T2 are upper triangular matrices whose eigenvalues lie */
	/*  along the diagonal.  The column indices ILO and IHI mark the starting */
	/*  and ending columns of the submatrix B. Balancing consists of applying */
	/*  a diagonal similarity transformation inv(D) * B * D to make the */
	/*  1-norms of each row of B and its corresponding column nearly equal. */
	/*  The output matrix is */

	/*     ( T1     X*D          Y    ) */
	/*     (  0  inv(D)*B*D  inv(D)*Z ). */
	/*     (  0      0           T2   ) */

	/*  Information about the permutations P and the diagonal matrix D is */
	/*  returned in the vector SCALE. */

	/*  This subroutine is based on the EISPACK routine BALANC. */

	/*  Modified by Tzu-Yi Chen, Computer Science Division, University of */
	/*    California at Berkeley, USA */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Test the input parameters */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--scale;

	/* Function Body */
	*info = 0;
	if (!lsame_ (job, "N", (ftnlen) 1, (ftnlen) 1) && !lsame_ (job, "P", (ftnlen) 1, (ftnlen) 1) && !lsame_ (job, "S", (ftnlen) 1, (ftnlen) 1)
		 && !lsame_ (job, "B", (ftnlen) 1, (ftnlen) 1))
	{
		*info = -1;
	}
	else if (*n < 0)
	{
		*info = -2;
	}
	else if (*lda < max (1, *n))
	{
		*info = -4;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DGEBAL", &i__1, (ftnlen) 6);
		return 0;
	}

	k = 1;
	l = *n;

	if (*n == 0)
	{
		goto L210;
	}

	if (lsame_ (job, "N", (ftnlen) 1, (ftnlen) 1))
	{
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			scale[i__] = 1.;
			/* L10: */
		}
		goto L210;
	}

	if (lsame_ (job, "S", (ftnlen) 1, (ftnlen) 1))
	{
		goto L120;
	}

	/*     Permutation to isolate eigenvalues if possible */

	goto L50;

	/*     Row and column exchange. */

 L20:
	scale[m] = (doublereal) j;
	if (j == m)
	{
		goto L30;
	}

	dswap_ (&l, &a[j * a_dim1 + 1], &c__1, &a[m * a_dim1 + 1], &c__1);
	i__1 = *n - k + 1;
	dswap_ (&i__1, &a[j + k * a_dim1], lda, &a[m + k * a_dim1], lda);

 L30:
	switch (iexc)
	{
	case 1:
		goto L40;
	case 2:
		goto L80;
	}

	/*     Search for rows isolating an eigenvalue and push them down. */

 L40:
	if (l == 1)
	{
		goto L210;
	}
	--l;

 L50:
	for (j = l; j >= 1; --j)
	{

		i__1 = l;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			if (i__ == j)
			{
				goto L60;
			}
			if (a[j + i__ * a_dim1] != 0.)
			{
				goto L70;
			}
		 L60:
			;
		}

		m = l;
		iexc = 1;
		goto L20;
	 L70:
		;
	}

	goto L90;

	/*     Search for columns isolating an eigenvalue and push them left. */

 L80:
	++k;

 L90:
	i__1 = l;
	for (j = k; j <= i__1; ++j)
	{

		i__2 = l;
		for (i__ = k; i__ <= i__2; ++i__)
		{
			if (i__ == j)
			{
				goto L100;
			}
			if (a[i__ + j * a_dim1] != 0.)
			{
				goto L110;
			}
		 L100:
			;
		}

		m = k;
		iexc = 2;
		goto L20;
	 L110:
		;
	}

 L120:
	i__1 = l;
	for (i__ = k; i__ <= i__1; ++i__)
	{
		scale[i__] = 1.;
		/* L130: */
	}

	if (lsame_ (job, "P", (ftnlen) 1, (ftnlen) 1))
	{
		goto L210;
	}

	/*     Balance the submatrix in rows K to L. */

	/*     Iterative loop for norm reduction */

	sfmin1 = dlamch_ ("S", (ftnlen) 1) / dlamch_ ("P", (ftnlen) 1);
	sfmax1 = 1. / sfmin1;
	sfmin2 = sfmin1 * 8.;
	sfmax2 = 1. / sfmin2;
 L140:
	noconv = FALSE_;

	i__1 = l;
	for (i__ = k; i__ <= i__1; ++i__)
	{
		c__ = 0.;
		r__ = 0.;

		i__2 = l;
		for (j = k; j <= i__2; ++j)
		{
			if (j == i__)
			{
				goto L150;
			}
			c__ += (d__1 = a[j + i__ * a_dim1], abs (d__1));
			r__ += (d__1 = a[i__ + j * a_dim1], abs (d__1));
		 L150:
			;
		}
		ica = idamax_ (&l, &a[i__ * a_dim1 + 1], &c__1);
		ca = (d__1 = a[ica + i__ * a_dim1], abs (d__1));
		i__2 = *n - k + 1;
		ira = idamax_ (&i__2, &a[i__ + k * a_dim1], lda);
		ra = (d__1 = a[i__ + (ira + k - 1) * a_dim1], abs (d__1));

		/*        Guard against zero C or R due to underflow. */

		if (c__ == 0. || r__ == 0.)
		{
			goto L200;
		}
		g = r__ / 8.;
		f = 1.;
		s = c__ + r__;
	 L160:
		/* Computing MAX */
		d__1 = max (f, c__);
		/* Computing MIN */
		d__2 = min (r__, g);
		if (c__ >= g || max (d__1, ca) >= sfmax2 || min (d__2, ra) <= sfmin2)
		{
			goto L170;
		}
		f *= 8.;
		c__ *= 8.;
		ca *= 8.;
		r__ /= 8.;
		g /= 8.;
		ra /= 8.;
		goto L160;

	 L170:
		g = c__ / 8.;
	 L180:
		/* Computing MIN */
		d__1 = min (f, c__), d__1 = min (d__1, g);
		if (g < r__ || max (r__, ra) >= sfmax2 || min (d__1, ca) <= sfmin2)
		{
			goto L190;
		}
		f /= 8.;
		c__ /= 8.;
		g /= 8.;
		ca /= 8.;
		r__ *= 8.;
		ra *= 8.;
		goto L180;

		/*        Now balance. */

	 L190:
		if (c__ + r__ >= s * .95)
		{
			goto L200;
		}
		if (f < 1. && scale[i__] < 1.)
		{
			if (f * scale[i__] <= sfmin1)
			{
				goto L200;
			}
		}
		if (f > 1. && scale[i__] > 1.)
		{
			if (scale[i__] >= sfmax1 / f)
			{
				goto L200;
			}
		}
		g = 1. / f;
		scale[i__] *= f;
		noconv = TRUE_;

		i__2 = *n - k + 1;
		dscal_ (&i__2, &g, &a[i__ + k * a_dim1], lda);
		dscal_ (&l, &f, &a[i__ * a_dim1 + 1], &c__1);

	 L200:
		;
	}

	if (noconv)
	{
		goto L140;
	}

 L210:
	*ilo = k;
	*ihi = l;

	return 0;

	/*     End of DGEBAL */

}	/* dgebal_ */

static int
dgecon_ (char *norm, integer * n, doublereal * a, integer *
			lda, doublereal * anorm, doublereal * rcond, doublereal * work, integer * iwork, integer * info, ftnlen norm_len)
{
	/* System generated locals */
	integer a_dim1, a_offset, i__1;
	doublereal d__1;

	/* Local variables */
	static doublereal sl;
	static integer ix;
	static doublereal su;
	static integer kase, kase1;
	static doublereal scale;

	static doublereal ainvnm;
	static logical onenrm;
	static char normin[1];
	static doublereal smlnum;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     February 29, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DGECON estimates the reciprocal of the condition number of a general */
	/*  real matrix A, in either the 1-norm or the infinity-norm, using */
	/*  the LU factorization computed by DGETRF. */

	/*  An estimate is obtained for norm(inv(A)), and the reciprocal of the */
	/*  condition number is computed as */
	/*     RCOND = 1 / ( norm(A) * norm(inv(A)) ). */

	/*  Arguments */
	/*  ========= */

	/*  NORM    (input) CHARACTER*1 */
	/*          Specifies whether the 1-norm condition number or the */
	/*          infinity-norm condition number is required: */
	/*          = '1' or 'O':  1-norm; */
	/*          = 'I':         Infinity-norm. */

	/*  N       (input) INTEGER */
	/*          The order of the matrix A.  N >= 0. */

	/*  A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          The factors L and U from the factorization A = P*L*U */
	/*          as computed by DGETRF. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(1,N). */

	/*  ANORM   (input) DOUBLE PRECISION */
	/*          If NORM = '1' or 'O', the 1-norm of the original matrix A. */
	/*          If NORM = 'I', the infinity-norm of the original matrix A. */

	/*  RCOND   (output) DOUBLE PRECISION */
	/*          The reciprocal of the condition number of the matrix A, */
	/*          computed as RCOND = 1/(norm(A) * norm(inv(A))). */

	/*  WORK    (workspace) DOUBLE PRECISION array, dimension (4*N) */

	/*  IWORK   (workspace) INTEGER array, dimension (N) */

	/*  INFO    (output) INTEGER */
	/*          = 0:  successful exit */
	/*          < 0:  if INFO = -i, the i-th argument had an illegal value */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Test the input parameters. */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--work;
	--iwork;

	/* Function Body */
	*info = 0;
	onenrm = *(unsigned char *) norm == '1' || lsame_ (norm, "O", (ftnlen) 1, (ftnlen) 1);
	if (!onenrm && !lsame_ (norm, "I", (ftnlen) 1, (ftnlen) 1))
	{
		*info = -1;
	}
	else if (*n < 0)
	{
		*info = -2;
	}
	else if (*lda < max (1, *n))
	{
		*info = -4;
	}
	else if (*anorm < 0.)
	{
		*info = -5;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DGECON", &i__1, (ftnlen) 6);
		return 0;
	}

	/*     Quick return if possible */

	*rcond = 0.;
	if (*n == 0)
	{
		*rcond = 1.;
		return 0;
	}
	else if (*anorm == 0.)
	{
		return 0;
	}

	smlnum = dlamch_ ("Safe minimum", (ftnlen) 12);

	/*     Estimate the norm of inv(A). */

	ainvnm = 0.;
	*(unsigned char *) normin = 'N';
	if (onenrm)
	{
		kase1 = 1;
	}
	else
	{
		kase1 = 2;
	}
	kase = 0;
 L10:
	dlacon_ (n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase);
	if (kase != 0)
	{
		if (kase == kase1)
		{

			/*           Multiply by inv(L). */

			dlatrs_ ("Lower", "No transpose", "Unit", normin, n, &a[a_offset],
						lda, &work[1], &sl, &work[(*n << 1) + 1], info, (ftnlen) 5, (ftnlen) 12, (ftnlen) 4, (ftnlen) 1);

			/*           Multiply by inv(U). */

			dlatrs_ ("Upper", "No transpose", "Non-unit", normin, n, &a[a_offset], lda, &work[1], &su, &work[*n * 3 + 1], info, (ftnlen) 5, (ftnlen) 12,
						(ftnlen) 8, (ftnlen) 1);
		}
		else
		{

			/*           Multiply by inv(U'). */

			dlatrs_ ("Upper", "Transpose", "Non-unit", normin, n, &a[a_offset],
						lda, &work[1], &su, &work[*n * 3 + 1], info, (ftnlen) 5, (ftnlen) 9, (ftnlen) 8, (ftnlen) 1);

			/*           Multiply by inv(L'). */

			dlatrs_ ("Lower", "Transpose", "Unit", normin, n, &a[a_offset],
						lda, &work[1], &sl, &work[(*n << 1) + 1], info, (ftnlen) 5, (ftnlen) 9, (ftnlen) 4, (ftnlen) 1);
		}

		/*        Divide X by 1/(SL*SU) if doing so will not cause overflow. */

		scale = sl * su;
		*(unsigned char *) normin = 'Y';
		if (scale != 1.)
		{
			ix = idamax_ (n, &work[1], &c__1);
			if (scale < (d__1 = work[ix], abs (d__1)) * smlnum || scale == 0.)
			{
				goto L20;
			}
			drscl_ (n, &scale, &work[1], &c__1);
		}
		goto L10;
	}

	/*     Compute the estimate of the reciprocal condition number. */

	if (ainvnm != 0.)
	{
		*rcond = 1. / ainvnm / *anorm;
	}

 L20:
	return 0;

	/*     End of DGECON */

}	/* dgecon_ */

static int
dgeequ_ (integer * m, integer * n, doublereal * a, integer *
			lda, doublereal * r__, doublereal * c__, doublereal * rowcnd, doublereal * colcnd, doublereal * amax, integer * info)
{
	/* System generated locals */
	integer a_dim1, a_offset, i__1, i__2;
	doublereal d__1, d__2, d__3;

	/* Local variables */
	static integer i__, j;
	static doublereal rcmin, rcmax;

	static doublereal bignum, smlnum;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     March 31, 1993 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DGEEQU computes row and column scalings intended to equilibrate an */
	/*  M-by-N matrix A and reduce its condition number.  R returns the row */
	/*  scale factors and C the column scale factors, chosen to try to make */
	/*  the largest element in each row and column of the matrix B with */
	/*  elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1. */

	/*  R(i) and C(j) are restricted to be between SMLNUM = smallest safe */
	/*  number and BIGNUM = largest safe number.  Use of these scaling */
	/*  factors is not guaranteed to reduce the condition number of A but */
	/*  works well in practice. */

	/*  Arguments */
	/*  ========= */

	/*  M       (input) INTEGER */
	/*          The number of rows of the matrix A.  M >= 0. */

	/*  N       (input) INTEGER */
	/*          The number of columns of the matrix A.  N >= 0. */

	/*  A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          The M-by-N matrix whose equilibration factors are */
	/*          to be computed. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(1,M). */

	/*  R       (output) DOUBLE PRECISION array, dimension (M) */
	/*          If INFO = 0 or INFO > M, R contains the row scale factors */
	/*          for A. */

	/*  C       (output) DOUBLE PRECISION array, dimension (N) */
	/*          If INFO = 0,  C contains the column scale factors for A. */

	/*  ROWCND  (output) DOUBLE PRECISION */
	/*          If INFO = 0 or INFO > M, ROWCND contains the ratio of the */
	/*          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and */
	/*          AMAX is neither too large nor too small, it is not worth */
	/*          scaling by R. */

	/*  COLCND  (output) DOUBLE PRECISION */
	/*          If INFO = 0, COLCND contains the ratio of the smallest */
	/*          C(i) to the largest C(i).  If COLCND >= 0.1, it is not */
	/*          worth scaling by C. */

	/*  AMAX    (output) DOUBLE PRECISION */
	/*          Absolute value of largest matrix element.  If AMAX is very */
	/*          close to overflow or very close to underflow, the matrix */
	/*          should be scaled. */

	/*  INFO    (output) INTEGER */
	/*          = 0:  successful exit */
	/*          < 0:  if INFO = -i, the i-th argument had an illegal value */
	/*          > 0:  if INFO = i,  and i is */
	/*                <= M:  the i-th row of A is exactly zero */
	/*                >  M:  the (i-M)-th column of A is exactly zero */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Test the input parameters. */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--r__;
	--c__;

	/* Function Body */
	*info = 0;
	if (*m < 0)
	{
		*info = -1;
	}
	else if (*n < 0)
	{
		*info = -2;
	}
	else if (*lda < max (1, *m))
	{
		*info = -4;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DGEEQU", &i__1, (ftnlen) 6);
		return 0;
	}

	/*     Quick return if possible */

	if (*m == 0 || *n == 0)
	{
		*rowcnd = 1.;
		*colcnd = 1.;
		*amax = 0.;
		return 0;
	}

	/*     Get machine constants. */

	smlnum = dlamch_ ("S", (ftnlen) 1);
	bignum = 1. / smlnum;

	/*     Compute row scale factors. */

	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		r__[i__] = 0.;
		/* L10: */
	}

	/*     Find the maximum element in each row. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j)
	{
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__)
		{
			/* Computing MAX */
			d__2 = r__[i__], d__3 = (d__1 = a[i__ + j * a_dim1], abs (d__1));
			r__[i__] = max (d__2, d__3);
			/* L20: */
		}
		/* L30: */
	}

	/*     Find the maximum and minimum scale factors. */

	rcmin = bignum;
	rcmax = 0.;
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		/* Computing MAX */
		d__1 = rcmax, d__2 = r__[i__];
		rcmax = max (d__1, d__2);
		/* Computing MIN */
		d__1 = rcmin, d__2 = r__[i__];
		rcmin = min (d__1, d__2);
		/* L40: */
	}
	*amax = rcmax;

	if (rcmin == 0.)
	{

		/*        Find the first zero scale factor and return an error code. */

		i__1 = *m;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			if (r__[i__] == 0.)
			{
				*info = i__;
				return 0;
			}
			/* L50: */
		}
	}
	else
	{

		/*        Invert the scale factors. */

		i__1 = *m;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			/* Computing MIN */
			/* Computing MAX */
			d__2 = r__[i__];
			d__1 = max (d__2, smlnum);
			r__[i__] = 1. / min (d__1, bignum);
			/* L60: */
		}

		/*        Compute ROWCND = min(R(I)) / max(R(I)) */

		*rowcnd = max (rcmin, smlnum) / min (rcmax, bignum);
	}

	/*     Compute column scale factors */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j)
	{
		c__[j] = 0.;
		/* L70: */
	}

	/*     Find the maximum element in each column, */
	/*     assuming the row scaling computed above. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j)
	{
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__)
		{
			/* Computing MAX */
			d__2 = c__[j], d__3 = (d__1 = a[i__ + j * a_dim1], abs (d__1)) * r__[i__];
			c__[j] = max (d__2, d__3);
			/* L80: */
		}
		/* L90: */
	}

	/*     Find the maximum and minimum scale factors. */

	rcmin = bignum;
	rcmax = 0.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j)
	{
		/* Computing MIN */
		d__1 = rcmin, d__2 = c__[j];
		rcmin = min (d__1, d__2);
		/* Computing MAX */
		d__1 = rcmax, d__2 = c__[j];
		rcmax = max (d__1, d__2);
		/* L100: */
	}

	if (rcmin == 0.)
	{

		/*        Find the first zero scale factor and return an error code. */

		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			if (c__[j] == 0.)
			{
				*info = *m + j;
				return 0;
			}
			/* L110: */
		}
	}
	else
	{

		/*        Invert the scale factors. */

		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			/* Computing MIN */
			/* Computing MAX */
			d__2 = c__[j];
			d__1 = max (d__2, smlnum);
			c__[j] = 1. / min (d__1, bignum);
			/* L120: */
		}

		/*        Compute COLCND = min(C(J)) / max(C(J)) */

		*colcnd = max (rcmin, smlnum) / min (rcmax, bignum);
	}

	return 0;

	/*     End of DGEEQU */

}	/* dgeequ_ */

int
dgeev_ (char *jobvl, char *jobvr, integer * n, doublereal *
		  a, integer * lda, doublereal * wr, doublereal * wi, doublereal * vl,
		  integer * ldvl, doublereal * vr, integer * ldvr, doublereal * work, integer * lwork, integer * info, ftnlen jobvl_len, ftnlen jobvr_len)
{
	/* System generated locals */
	integer a_dim1, a_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, i__2, i__3, i__4;
	doublereal d__1, d__2;

	/* Builtin functions */
	double sqrt (doublereal);

	/* Local variables */
	static integer i__, k;
	static doublereal r__, cs, sn;
	static integer ihi;
	static doublereal scl;
	static integer ilo;
	static doublereal dum[1], eps;
	static integer ibal;
	static char side[1];
	static integer maxb;
	static doublereal anrm;
	static integer ierr, itau;

	static integer iwrk, nout;

	static logical scalea;
	static doublereal cscale;
	static logical select[1];

	static doublereal bignum;

	static integer minwrk, maxwrk;
	static logical wantvl;
	static doublereal smlnum;
	static integer hswork;
	static logical lquery, wantvr;


	/*  -- LAPACK driver routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     December 8, 1999 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DGEEV computes for an N-by-N real nonsymmetric matrix A, the */
	/*  eigenvalues and, optionally, the left and/or right eigenvectors. */

	/*  The right eigenvector v(j) of A satisfies */
	/*                   A * v(j) = lambda(j) * v(j) */
	/*  where lambda(j) is its eigenvalue. */
	/*  The left eigenvector u(j) of A satisfies */
	/*                u(j)**H * A = lambda(j) * u(j)**H */
	/*  where u(j)**H denotes the conjugate transpose of u(j). */

	/*  The computed eigenvectors are normalized to have Euclidean norm */
	/*  equal to 1 and largest component real. */

	/*  Arguments */
	/*  ========= */

	/*  JOBVL   (input) CHARACTER*1 */
	/*          = 'N': left eigenvectors of A are not computed; */
	/*          = 'V': left eigenvectors of A are computed. */

	/*  JOBVR   (input) CHARACTER*1 */
	/*          = 'N': right eigenvectors of A are not computed; */
	/*          = 'V': right eigenvectors of A are computed. */

	/*  N       (input) INTEGER */
	/*          The order of the matrix A. N >= 0. */

	/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          On entry, the N-by-N matrix A. */
	/*          On exit, A has been overwritten. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(1,N). */

	/*  WR      (output) DOUBLE PRECISION array, dimension (N) */
	/*  WI      (output) DOUBLE PRECISION array, dimension (N) */
	/*          WR and WI contain the real and imaginary parts, */
	/*          respectively, of the computed eigenvalues.  Complex */
	/*          conjugate pairs of eigenvalues appear consecutively */
	/*          with the eigenvalue having the positive imaginary part */
	/*          first. */

	/*  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N) */
	/*          If JOBVL = 'V', the left eigenvectors u(j) are stored one */
	/*          after another in the columns of VL, in the same order */
	/*          as their eigenvalues. */
	/*          If JOBVL = 'N', VL is not referenced. */
	/*          If the j-th eigenvalue is real, then u(j) = VL(:,j), */
	/*          the j-th column of VL. */
	/*          If the j-th and (j+1)-st eigenvalues form a complex */
	/*          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and */
	/*          u(j+1) = VL(:,j) - i*VL(:,j+1). */

	/*  LDVL    (input) INTEGER */
	/*          The leading dimension of the array VL.  LDVL >= 1; if */
	/*          JOBVL = 'V', LDVL >= N. */

	/*  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N) */
	/*          If JOBVR = 'V', the right eigenvectors v(j) are stored one */
	/*          after another in the columns of VR, in the same order */
	/*          as their eigenvalues. */
	/*          If JOBVR = 'N', VR is not referenced. */
	/*          If the j-th eigenvalue is real, then v(j) = VR(:,j), */
	/*          the j-th column of VR. */
	/*          If the j-th and (j+1)-st eigenvalues form a complex */
	/*          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and */
	/*          v(j+1) = VR(:,j) - i*VR(:,j+1). */

	/*  LDVR    (input) INTEGER */
	/*          The leading dimension of the array VR.  LDVR >= 1; if */
	/*          JOBVR = 'V', LDVR >= N. */

	/*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) */
	/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

	/*  LWORK   (input) INTEGER */
	/*          The dimension of the array WORK.  LWORK >= max(1,3*N), and */
	/*          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good */
	/*          performance, LWORK must generally be larger. */

	/*          If LWORK = -1, then a workspace query is assumed; the routine */
	/*          only calculates the optimal size of the WORK array, returns */
	/*          this value as the first entry of the WORK array, and no error */
	/*          message related to LWORK is issued by XERBLA. */

	/*  INFO    (output) INTEGER */
	/*          = 0:  successful exit */
	/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */
	/*          > 0:  if INFO = i, the QR algorithm failed to compute all the */
	/*                eigenvalues, and no eigenvectors have been computed; */
	/*                elements i+1:N of WR and WI contain eigenvalues which */
	/*                have converged. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Local Arrays .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Test the input arguments */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--wr;
	--wi;
	vl_dim1 = *ldvl;
	vl_offset = 1 + vl_dim1;
	vl -= vl_offset;
	vr_dim1 = *ldvr;
	vr_offset = 1 + vr_dim1;
	vr -= vr_offset;
	--work;

	/* Function Body */
	*info = 0;
	lquery = *lwork == -1;
	wantvl = lsame_ (jobvl, "V", (ftnlen) 1, (ftnlen) 1);
	wantvr = lsame_ (jobvr, "V", (ftnlen) 1, (ftnlen) 1);
	if (!wantvl && !lsame_ (jobvl, "N", (ftnlen) 1, (ftnlen) 1))
	{
		*info = -1;
	}
	else if (!wantvr && !lsame_ (jobvr, "N", (ftnlen) 1, (ftnlen) 1))
	{
		*info = -2;
	}
	else if (*n < 0)
	{
		*info = -3;
	}
	else if (*lda < max (1, *n))
	{
		*info = -5;
	}
	else if (*ldvl < 1 || wantvl && *ldvl < *n)
	{
		*info = -9;
	}
	else if (*ldvr < 1 || wantvr && *ldvr < *n)
	{
		*info = -11;
	}

	/*     Compute workspace */
	/*      (Note: Comments in the code beginning "Workspace:" describe the */
	/*       minimal amount of workspace needed at that point in the code, */
	/*       as well as the preferred amount for good performance. */
	/*       NB refers to the optimal block size for the immediately */
	/*       following subroutine, as returned by ILAENV. */
	/*       HSWORK refers to the workspace preferred by DHSEQR, as */
	/*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N, */
	/*       the worst case.) */

	minwrk = 1;
	if (*info == 0 && (*lwork >= 1 || lquery))
	{
		maxwrk = (*n << 1) + *n * ilaenv_ (&c__1, "DGEHRD", " ", n, &c__1, n, &c__0, (ftnlen) 6, (ftnlen) 1);
		if (!wantvl && !wantvr)
		{
			/* Computing MAX */
			i__1 = 1, i__2 = *n * 3;
			minwrk = max (i__1, i__2);
			/* Computing MAX */
			i__1 = ilaenv_ (&c__8, "DHSEQR", "EN", n, &c__1, n, &c_n1, (ftnlen) 6, (ftnlen) 2);
			maxb = max (i__1, 2);
			/* Computing MIN */
			/* Computing MAX */
			i__3 = 2, i__4 = ilaenv_ (&c__4, "DHSEQR", "EN", n, &c__1, n, &c_n1, (ftnlen) 6, (ftnlen) 2);
			i__1 = min (maxb, *n), i__2 = max (i__3, i__4);
			k = min (i__1, i__2);
			/* Computing MAX */
			i__1 = k * (k + 2), i__2 = *n << 1;
			hswork = max (i__1, i__2);
			/* Computing MAX */
			i__1 = maxwrk, i__2 = *n + 1, i__1 = max (i__1, i__2), i__2 = *n + hswork;
			maxwrk = max (i__1, i__2);
		}
		else
		{
			/* Computing MAX */
			i__1 = 1, i__2 = *n << 2;
			minwrk = max (i__1, i__2);
			/* Computing MAX */
			i__1 = maxwrk, i__2 = (*n << 1) + (*n - 1) * ilaenv_ (&c__1, "DOR" "GHR", " ", n, &c__1, n, &c_n1, (ftnlen) 6, (ftnlen) 1);
			maxwrk = max (i__1, i__2);
			/* Computing MAX */
			i__1 = ilaenv_ (&c__8, "DHSEQR", "SV", n, &c__1, n, &c_n1, (ftnlen) 6, (ftnlen) 2);
			maxb = max (i__1, 2);
			/* Computing MIN */
			/* Computing MAX */
			i__3 = 2, i__4 = ilaenv_ (&c__4, "DHSEQR", "SV", n, &c__1, n, &c_n1, (ftnlen) 6, (ftnlen) 2);
			i__1 = min (maxb, *n), i__2 = max (i__3, i__4);
			k = min (i__1, i__2);
			/* Computing MAX */
			i__1 = k * (k + 2), i__2 = *n << 1;
			hswork = max (i__1, i__2);
			/* Computing MAX */
			i__1 = maxwrk, i__2 = *n + 1, i__1 = max (i__1, i__2), i__2 = *n + hswork;
			maxwrk = max (i__1, i__2);
			/* Computing MAX */
			i__1 = maxwrk, i__2 = *n << 2;
			maxwrk = max (i__1, i__2);
		}
		work[1] = (doublereal) maxwrk;
	}
	if (*lwork < minwrk && !lquery)
	{
		*info = -13;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DGEEV ", &i__1, (ftnlen) 6);
		return 0;
	}
	else if (lquery)
	{
		return 0;
	}

	/*     Quick return if possible */

	if (*n == 0)
	{
		return 0;
	}

	/*     Get machine constants */

	eps = dlamch_ ("P", (ftnlen) 1);
	smlnum = dlamch_ ("S", (ftnlen) 1);
	bignum = 1. / smlnum;
	dlabad_ (&smlnum, &bignum);
	smlnum = sqrt (smlnum) / eps;
	bignum = 1. / smlnum;

	/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

	anrm = dlange_ ("M", n, n, &a[a_offset], lda, dum, (ftnlen) 1);
	scalea = FALSE_;
	if (anrm > 0. && anrm < smlnum)
	{
		scalea = TRUE_;
		cscale = smlnum;
	}
	else if (anrm > bignum)
	{
		scalea = TRUE_;
		cscale = bignum;
	}
	if (scalea)
	{
		dlascl_ ("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &ierr, (ftnlen) 1);
	}

	/*     Balance the matrix */
	/*     (Workspace: need N) */

	ibal = 1;
	dgebal_ ("B", n, &a[a_offset], lda, &ilo, &ihi, &work[ibal], &ierr, (ftnlen) 1);

	/*     Reduce to upper Hessenberg form */
	/*     (Workspace: need 3*N, prefer 2*N+N*NB) */

	itau = ibal + *n;
	iwrk = itau + *n;
	i__1 = *lwork - iwrk + 1;
	dgehrd_ (n, &ilo, &ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1, &ierr);

	if (wantvl)
	{

		/*        Want left eigenvectors */
		/*        Copy Householder vectors to VL */

		*(unsigned char *) side = 'L';
		dlacpy_ ("L", n, n, &a[a_offset], lda, &vl[vl_offset], ldvl, (ftnlen) 1);

		/*        Generate orthogonal matrix in VL */
		/*        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB) */

		i__1 = *lwork - iwrk + 1;
		dorghr_ (n, &ilo, &ihi, &vl[vl_offset], ldvl, &work[itau], &work[iwrk], &i__1, &ierr);

		/*        Perform QR iteration, accumulating Schur vectors in VL */
		/*        (Workspace: need N+1, prefer N+HSWORK (see comments) ) */

		iwrk = itau;
		i__1 = *lwork - iwrk + 1;
		dhseqr_ ("S", "V", n, &ilo, &ihi, &a[a_offset], lda, &wr[1], &wi[1], &vl[vl_offset], ldvl, &work[iwrk], &i__1, info, (ftnlen) 1, (ftnlen) 1);

		if (wantvr)
		{

			/*           Want left and right eigenvectors */
			/*           Copy Schur vectors to VR */

			*(unsigned char *) side = 'B';
			dlacpy_ ("F", n, n, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, (ftnlen) 1);
		}

	}
	else if (wantvr)
	{

		/*        Want right eigenvectors */
		/*        Copy Householder vectors to VR */

		*(unsigned char *) side = 'R';
		dlacpy_ ("L", n, n, &a[a_offset], lda, &vr[vr_offset], ldvr, (ftnlen) 1);

		/*        Generate orthogonal matrix in VR */
		/*        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB) */

		i__1 = *lwork - iwrk + 1;
		dorghr_ (n, &ilo, &ihi, &vr[vr_offset], ldvr, &work[itau], &work[iwrk], &i__1, &ierr);

		/*        Perform QR iteration, accumulating Schur vectors in VR */
		/*        (Workspace: need N+1, prefer N+HSWORK (see comments) ) */

		iwrk = itau;
		i__1 = *lwork - iwrk + 1;
		dhseqr_ ("S", "V", n, &ilo, &ihi, &a[a_offset], lda, &wr[1], &wi[1], &vr[vr_offset], ldvr, &work[iwrk], &i__1, info, (ftnlen) 1, (ftnlen) 1);

	}
	else
	{

		/*        Compute eigenvalues only */
		/*        (Workspace: need N+1, prefer N+HSWORK (see comments) ) */

		iwrk = itau;
		i__1 = *lwork - iwrk + 1;
		dhseqr_ ("E", "N", n, &ilo, &ihi, &a[a_offset], lda, &wr[1], &wi[1], &vr[vr_offset], ldvr, &work[iwrk], &i__1, info, (ftnlen) 1, (ftnlen) 1);
	}

	/*     If INFO > 0 from DHSEQR, then quit */

	if (*info > 0)
	{
		goto L50;
	}

	if (wantvl || wantvr)
	{

		/*        Compute left and/or right eigenvectors */
		/*        (Workspace: need 4*N) */

		dtrevc_ (side, "B", select, n, &a[a_offset], lda, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, n, &nout, &work[iwrk], &ierr, (ftnlen) 1, (ftnlen) 1);
	}

	if (wantvl)
	{

		/*        Undo balancing of left eigenvectors */
		/*        (Workspace: need N) */

		dgebak_ ("B", "L", n, &ilo, &ihi, &work[ibal], n, &vl[vl_offset], ldvl, &ierr, (ftnlen) 1, (ftnlen) 1);

		/*        Normalize left eigenvectors and make largest component real */

		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			if (wi[i__] == 0.)
			{
				scl = 1. / dnrm2_ (n, &vl[i__ * vl_dim1 + 1], &c__1);
				dscal_ (n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
			}
			else if (wi[i__] > 0.)
			{
				d__1 = dnrm2_ (n, &vl[i__ * vl_dim1 + 1], &c__1);
				d__2 = dnrm2_ (n, &vl[(i__ + 1) * vl_dim1 + 1], &c__1);
				scl = 1. / dlapy2_ (&d__1, &d__2);
				dscal_ (n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
				dscal_ (n, &scl, &vl[(i__ + 1) * vl_dim1 + 1], &c__1);
				i__2 = *n;
				for (k = 1; k <= i__2; ++k)
				{
					/* Computing 2nd power */
					d__1 = vl[k + i__ * vl_dim1];
					/* Computing 2nd power */
					d__2 = vl[k + (i__ + 1) * vl_dim1];
					work[iwrk + k - 1] = d__1 * d__1 + d__2 * d__2;
					/* L10: */
				}
				k = idamax_ (n, &work[iwrk], &c__1);
				dlartg_ (&vl[k + i__ * vl_dim1], &vl[k + (i__ + 1) * vl_dim1], &cs, &sn, &r__);
				drot_ (n, &vl[i__ * vl_dim1 + 1], &c__1, &vl[(i__ + 1) * vl_dim1 + 1], &c__1, &cs, &sn);
				vl[k + (i__ + 1) * vl_dim1] = 0.;
			}
			/* L20: */
		}
	}

	if (wantvr)
	{

		/*        Undo balancing of right eigenvectors */
		/*        (Workspace: need N) */

		dgebak_ ("B", "R", n, &ilo, &ihi, &work[ibal], n, &vr[vr_offset], ldvr, &ierr, (ftnlen) 1, (ftnlen) 1);

		/*        Normalize right eigenvectors and make largest component real */

		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			if (wi[i__] == 0.)
			{
				scl = 1. / dnrm2_ (n, &vr[i__ * vr_dim1 + 1], &c__1);
				dscal_ (n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
			}
			else if (wi[i__] > 0.)
			{
				d__1 = dnrm2_ (n, &vr[i__ * vr_dim1 + 1], &c__1);
				d__2 = dnrm2_ (n, &vr[(i__ + 1) * vr_dim1 + 1], &c__1);
				scl = 1. / dlapy2_ (&d__1, &d__2);
				dscal_ (n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
				dscal_ (n, &scl, &vr[(i__ + 1) * vr_dim1 + 1], &c__1);
				i__2 = *n;
				for (k = 1; k <= i__2; ++k)
				{
					/* Computing 2nd power */
					d__1 = vr[k + i__ * vr_dim1];
					/* Computing 2nd power */
					d__2 = vr[k + (i__ + 1) * vr_dim1];
					work[iwrk + k - 1] = d__1 * d__1 + d__2 * d__2;
					/* L30: */
				}
				k = idamax_ (n, &work[iwrk], &c__1);
				dlartg_ (&vr[k + i__ * vr_dim1], &vr[k + (i__ + 1) * vr_dim1], &cs, &sn, &r__);
				drot_ (n, &vr[i__ * vr_dim1 + 1], &c__1, &vr[(i__ + 1) * vr_dim1 + 1], &c__1, &cs, &sn);
				vr[k + (i__ + 1) * vr_dim1] = 0.;
			}
			/* L40: */
		}
	}

	/*     Undo scaling if necessary */

 L50:
	if (scalea)
	{
		i__1 = *n - *info;
		/* Computing MAX */
		i__3 = *n - *info;
		i__2 = max (i__3, 1);
		dlascl_ ("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wr[*info + 1], &i__2, &ierr, (ftnlen) 1);
		i__1 = *n - *info;
		/* Computing MAX */
		i__3 = *n - *info;
		i__2 = max (i__3, 1);
		dlascl_ ("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[*info + 1], &i__2, &ierr, (ftnlen) 1);
		if (*info > 0)
		{
			i__1 = ilo - 1;
			dlascl_ ("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wr[1], n, &ierr, (ftnlen) 1);
			i__1 = ilo - 1;
			dlascl_ ("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[1], n, &ierr, (ftnlen) 1);
		}
	}

	work[1] = (doublereal) maxwrk;
	return 0;

	/*     End of DGEEV */

}	/* dgeev_ */

int
dgeevx_ (char *balanc, char *jobvl, char *jobvr, char *sense, integer * n, doublereal * a, integer * lda, doublereal * wr,
			doublereal * wi, doublereal * vl, integer * ldvl, doublereal * vr,
			integer * ldvr, integer * ilo, integer * ihi, doublereal * scale,
			doublereal * abnrm, doublereal * rconde, doublereal * rcondv, doublereal
			* work, integer * lwork, integer * iwork, integer * info, ftnlen balanc_len, ftnlen jobvl_len, ftnlen jobvr_len, ftnlen sense_len)
{
	/* System generated locals */
	integer a_dim1, a_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, i__2, i__3, i__4;
	doublereal d__1, d__2;

	/* Builtin functions */
	double sqrt (doublereal);

	/* Local variables */
	static integer i__, k;
	static doublereal r__, cs, sn;
	static char job[1];
	static doublereal scl, dum[1], eps;
	static char side[1];
	static integer maxb;
	static doublereal anrm;
	static integer ierr, itau;

	static integer iwrk, nout;
	static integer icond;
	static logical scalea;
	static doublereal cscale;
	static logical select[1];
	static doublereal bignum;
	static integer minwrk, maxwrk;
	static logical wantvl, wntsnb;
	static integer hswork;
	static logical wntsne;
	static doublereal smlnum;
	static logical lquery, wantvr, wntsnn, wntsnv;


	/*  -- LAPACK driver routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     June 30, 1999 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DGEEVX computes for an N-by-N real nonsymmetric matrix A, the */
	/*  eigenvalues and, optionally, the left and/or right eigenvectors. */

	/*  Optionally also, it computes a balancing transformation to improve */
	/*  the conditioning of the eigenvalues and eigenvectors (ILO, IHI, */
	/*  SCALE, and ABNRM), reciprocal condition numbers for the eigenvalues */
	/*  (RCONDE), and reciprocal condition numbers for the right */
	/*  eigenvectors (RCONDV). */

	/*  The right eigenvector v(j) of A satisfies */
	/*                   A * v(j) = lambda(j) * v(j) */
	/*  where lambda(j) is its eigenvalue. */
	/*  The left eigenvector u(j) of A satisfies */
	/*                u(j)**H * A = lambda(j) * u(j)**H */
	/*  where u(j)**H denotes the conjugate transpose of u(j). */

	/*  The computed eigenvectors are normalized to have Euclidean norm */
	/*  equal to 1 and largest component real. */

	/*  Balancing a matrix means permuting the rows and columns to make it */
	/*  more nearly upper triangular, and applying a diagonal similarity */
	/*  transformation D * A * D**(-1), where D is a diagonal matrix, to */
	/*  make its rows and columns closer in norm and the condition numbers */
	/*  of its eigenvalues and eigenvectors smaller.  The computed */
	/*  reciprocal condition numbers correspond to the balanced matrix. */
	/*  Permuting rows and columns will not change the condition numbers */
	/*  (in exact arithmetic) but diagonal scaling will.  For further */
	/*  explanation of balancing, see section 4.10.2 of the LAPACK */
	/*  Users' Guide. */

	/*  Arguments */
	/*  ========= */

	/*  BALANC  (input) CHARACTER*1 */
	/*          Indicates how the input matrix should be diagonally scaled */
	/*          and/or permuted to improve the conditioning of its */
	/*          eigenvalues. */
	/*          = 'N': Do not diagonally scale or permute; */
	/*          = 'P': Perform permutations to make the matrix more nearly */
	/*                 upper triangular. Do not diagonally scale; */
	/*          = 'S': Diagonally scale the matrix, i.e. replace A by */
	/*                 D*A*D**(-1), where D is a diagonal matrix chosen */
	/*                 to make the rows and columns of A more equal in */
	/*                 norm. Do not permute; */
	/*          = 'B': Both diagonally scale and permute A. */

	/*          Computed reciprocal condition numbers will be for the matrix */
	/*          after balancing and/or permuting. Permuting does not change */
	/*          condition numbers (in exact arithmetic), but balancing does. */

	/*  JOBVL   (input) CHARACTER*1 */
	/*          = 'N': left eigenvectors of A are not computed; */
	/*          = 'V': left eigenvectors of A are computed. */
	/*          If SENSE = 'E' or 'B', JOBVL must = 'V'. */

	/*  JOBVR   (input) CHARACTER*1 */
	/*          = 'N': right eigenvectors of A are not computed; */
	/*          = 'V': right eigenvectors of A are computed. */
	/*          If SENSE = 'E' or 'B', JOBVR must = 'V'. */

	/*  SENSE   (input) CHARACTER*1 */
	/*          Determines which reciprocal condition numbers are computed. */
	/*          = 'N': None are computed; */
	/*          = 'E': Computed for eigenvalues only; */
	/*          = 'V': Computed for right eigenvectors only; */
	/*          = 'B': Computed for eigenvalues and right eigenvectors. */

	/*          If SENSE = 'E' or 'B', both left and right eigenvectors */
	/*          must also be computed (JOBVL = 'V' and JOBVR = 'V'). */

	/*  N       (input) INTEGER */
	/*          The order of the matrix A. N >= 0. */

	/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          On entry, the N-by-N matrix A. */
	/*          On exit, A has been overwritten.  If JOBVL = 'V' or */
	/*          JOBVR = 'V', A contains the real Schur form of the balanced */
	/*          version of the input matrix A. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(1,N). */

	/*  WR      (output) DOUBLE PRECISION array, dimension (N) */
	/*  WI      (output) DOUBLE PRECISION array, dimension (N) */
	/*          WR and WI contain the real and imaginary parts, */
	/*          respectively, of the computed eigenvalues.  Complex */
	/*          conjugate pairs of eigenvalues will appear consecutively */
	/*          with the eigenvalue having the positive imaginary part */
	/*          first. */

	/*  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N) */
	/*          If JOBVL = 'V', the left eigenvectors u(j) are stored one */
	/*          after another in the columns of VL, in the same order */
	/*          as their eigenvalues. */
	/*          If JOBVL = 'N', VL is not referenced. */
	/*          If the j-th eigenvalue is real, then u(j) = VL(:,j), */
	/*          the j-th column of VL. */
	/*          If the j-th and (j+1)-st eigenvalues form a complex */
	/*          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and */
	/*          u(j+1) = VL(:,j) - i*VL(:,j+1). */

	/*  LDVL    (input) INTEGER */
	/*          The leading dimension of the array VL.  LDVL >= 1; if */
	/*          JOBVL = 'V', LDVL >= N. */

	/*  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N) */
	/*          If JOBVR = 'V', the right eigenvectors v(j) are stored one */
	/*          after another in the columns of VR, in the same order */
	/*          as their eigenvalues. */
	/*          If JOBVR = 'N', VR is not referenced. */
	/*          If the j-th eigenvalue is real, then v(j) = VR(:,j), */
	/*          the j-th column of VR. */
	/*          If the j-th and (j+1)-st eigenvalues form a complex */
	/*          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and */
	/*          v(j+1) = VR(:,j) - i*VR(:,j+1). */

	/*  LDVR    (input) INTEGER */
	/*          The leading dimension of the array VR.  LDVR >= 1, and if */
	/*          JOBVR = 'V', LDVR >= N. */

	/*  ILO,IHI (output) INTEGER */
	/*          ILO and IHI are integer values determined when A was */
	/*          balanced.  The balanced A(i,j) = 0 if I > J and */
	/*          J = 1,...,ILO-1 or I = IHI+1,...,N. */

	/*  SCALE   (output) DOUBLE PRECISION array, dimension (N) */
	/*          Details of the permutations and scaling factors applied */
	/*          when balancing A.  If P(j) is the index of the row and column */
	/*          interchanged with row and column j, and D(j) is the scaling */
	/*          factor applied to row and column j, then */
	/*          SCALE(J) = P(J),    for J = 1,...,ILO-1 */
	/*                   = D(J),    for J = ILO,...,IHI */
	/*                   = P(J)     for J = IHI+1,...,N. */
	/*          The order in which the interchanges are made is N to IHI+1, */
	/*          then 1 to ILO-1. */

	/*  ABNRM   (output) DOUBLE PRECISION */
	/*          The one-norm of the balanced matrix (the maximum */
	/*          of the sum of absolute values of elements of any column). */

	/*  RCONDE  (output) DOUBLE PRECISION array, dimension (N) */
	/*          RCONDE(j) is the reciprocal condition number of the j-th */
	/*          eigenvalue. */

	/*  RCONDV  (output) DOUBLE PRECISION array, dimension (N) */
	/*          RCONDV(j) is the reciprocal condition number of the j-th */
	/*          right eigenvector. */

	/*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) */
	/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

	/*  LWORK   (input) INTEGER */
	/*          The dimension of the array WORK.   If SENSE = 'N' or 'E', */
	/*          LWORK >= max(1,2*N), and if JOBVL = 'V' or JOBVR = 'V', */
	/*          LWORK >= 3*N.  If SENSE = 'V' or 'B', LWORK >= N*(N+6). */
	/*          For good performance, LWORK must generally be larger. */

	/*          If LWORK = -1, then a workspace query is assumed; the routine */
	/*          only calculates the optimal size of the WORK array, returns */
	/*          this value as the first entry of the WORK array, and no error */
	/*          message related to LWORK is issued by XERBLA. */

	/*  IWORK   (workspace) INTEGER array, dimension (2*N-2) */
	/*          If SENSE = 'N' or 'E', not referenced. */

	/*  INFO    (output) INTEGER */
	/*          = 0:  successful exit */
	/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */
	/*          > 0:  if INFO = i, the QR algorithm failed to compute all the */
	/*                eigenvalues, and no eigenvectors or condition numbers */
	/*                have been computed; elements 1:ILO-1 and i+1:N of WR */
	/*                and WI contain eigenvalues which have converged. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Local Arrays .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Test the input arguments */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--wr;
	--wi;
	vl_dim1 = *ldvl;
	vl_offset = 1 + vl_dim1;
	vl -= vl_offset;
	vr_dim1 = *ldvr;
	vr_offset = 1 + vr_dim1;
	vr -= vr_offset;
	--scale;
	--rconde;
	--rcondv;
	--work;
	--iwork;

	/* Function Body */
	*info = 0;
	lquery = *lwork == -1;
	wantvl = lsame_ (jobvl, "V", (ftnlen) 1, (ftnlen) 1);
	wantvr = lsame_ (jobvr, "V", (ftnlen) 1, (ftnlen) 1);
	wntsnn = lsame_ (sense, "N", (ftnlen) 1, (ftnlen) 1);
	wntsne = lsame_ (sense, "E", (ftnlen) 1, (ftnlen) 1);
	wntsnv = lsame_ (sense, "V", (ftnlen) 1, (ftnlen) 1);
	wntsnb = lsame_ (sense, "B", (ftnlen) 1, (ftnlen) 1);
	if (!(lsame_ (balanc, "N", (ftnlen) 1, (ftnlen) 1) || lsame_ (balanc, "S", (ftnlen) 1, (ftnlen) 1) || lsame_ (balanc, "P", (ftnlen) 1, (ftnlen) 1)
			|| lsame_ (balanc, "B", (ftnlen) 1, (ftnlen) 1)))
	{
		*info = -1;
	}
	else if (!wantvl && !lsame_ (jobvl, "N", (ftnlen) 1, (ftnlen) 1))
	{
		*info = -2;
	}
	else if (!wantvr && !lsame_ (jobvr, "N", (ftnlen) 1, (ftnlen) 1))
	{
		*info = -3;
	}
	else if (!(wntsnn || wntsne || wntsnb || wntsnv) || (wntsne || wntsnb) && !(wantvl && wantvr))
	{
		*info = -4;
	}
	else if (*n < 0)
	{
		*info = -5;
	}
	else if (*lda < max (1, *n))
	{
		*info = -7;
	}
	else if (*ldvl < 1 || wantvl && *ldvl < *n)
	{
		*info = -11;
	}
	else if (*ldvr < 1 || wantvr && *ldvr < *n)
	{
		*info = -13;
	}

	/*     Compute workspace */
	/*      (Note: Comments in the code beginning "Workspace:" describe the */
	/*       minimal amount of workspace needed at that point in the code, */
	/*       as well as the preferred amount for good performance. */
	/*       NB refers to the optimal block size for the immediately */
	/*       following subroutine, as returned by ILAENV. */
	/*       HSWORK refers to the workspace preferred by DHSEQR, as */
	/*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N, */
	/*       the worst case.) */

	minwrk = 1;
	if (*info == 0 && (*lwork >= 1 || lquery))
	{
		maxwrk = *n + *n * ilaenv_ (&c__1, "DGEHRD", " ", n, &c__1, n, &c__0, (ftnlen) 6, (ftnlen) 1);
		if (!wantvl && !wantvr)
		{
			/* Computing MAX */
			i__1 = 1, i__2 = *n << 1;
			minwrk = max (i__1, i__2);
			if (!wntsnn)
			{
				/* Computing MAX */
				i__1 = minwrk, i__2 = *n * *n + *n * 6;
				minwrk = max (i__1, i__2);
			}
			/* Computing MAX */
			i__1 = ilaenv_ (&c__8, "DHSEQR", "SN", n, &c__1, n, &c_n1, (ftnlen) 6, (ftnlen) 2);
			maxb = max (i__1, 2);
			if (wntsnn)
			{
				/* Computing MIN */
				/* Computing MAX */
				i__3 = 2, i__4 = ilaenv_ (&c__4, "DHSEQR", "EN", n, &c__1, n, &c_n1, (ftnlen) 6, (ftnlen) 2);
				i__1 = min (maxb, *n), i__2 = max (i__3, i__4);
				k = min (i__1, i__2);
			}
			else
			{
				/* Computing MIN */
				/* Computing MAX */
				i__3 = 2, i__4 = ilaenv_ (&c__4, "DHSEQR", "SN", n, &c__1, n, &c_n1, (ftnlen) 6, (ftnlen) 2);
				i__1 = min (maxb, *n), i__2 = max (i__3, i__4);
				k = min (i__1, i__2);
			}
			/* Computing MAX */
			i__1 = k * (k + 2), i__2 = *n << 1;
			hswork = max (i__1, i__2);
			/* Computing MAX */
			i__1 = max (maxwrk, 1);
			maxwrk = max (i__1, hswork);
			if (!wntsnn)
			{
				/* Computing MAX */
				i__1 = maxwrk, i__2 = *n * *n + *n * 6;
				maxwrk = max (i__1, i__2);
			}
		}
		else
		{
			/* Computing MAX */
			i__1 = 1, i__2 = *n * 3;
			minwrk = max (i__1, i__2);
			if (!wntsnn && !wntsne)
			{
				/* Computing MAX */
				i__1 = minwrk, i__2 = *n * *n + *n * 6;
				minwrk = max (i__1, i__2);
			}
			/* Computing MAX */
			i__1 = ilaenv_ (&c__8, "DHSEQR", "SN", n, &c__1, n, &c_n1, (ftnlen) 6, (ftnlen) 2);
			maxb = max (i__1, 2);
			/* Computing MIN */
			/* Computing MAX */
			i__3 = 2, i__4 = ilaenv_ (&c__4, "DHSEQR", "EN", n, &c__1, n, &c_n1, (ftnlen) 6, (ftnlen) 2);
			i__1 = min (maxb, *n), i__2 = max (i__3, i__4);
			k = min (i__1, i__2);
			/* Computing MAX */
			i__1 = k * (k + 2), i__2 = *n << 1;
			hswork = max (i__1, i__2);
			/* Computing MAX */
			i__1 = max (maxwrk, 1);
			maxwrk = max (i__1, hswork);
			/* Computing MAX */
			i__1 = maxwrk, i__2 = *n + (*n - 1) * ilaenv_ (&c__1, "DORGHR", " ", n, &c__1, n, &c_n1, (ftnlen) 6, (ftnlen) 1);
			maxwrk = max (i__1, i__2);
			if (!wntsnn && !wntsne)
			{
				/* Computing MAX */
				i__1 = maxwrk, i__2 = *n * *n + *n * 6;
				maxwrk = max (i__1, i__2);
			}
			/* Computing MAX */
			i__1 = maxwrk, i__2 = *n * 3, i__1 = max (i__1, i__2);
			maxwrk = max (i__1, 1);
		}
		work[1] = (doublereal) maxwrk;
	}
	if (*lwork < minwrk && !lquery)
	{
		*info = -21;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DGEEVX", &i__1, (ftnlen) 6);
		return 0;
	}
	else if (lquery)
	{
		return 0;
	}

	/*     Quick return if possible */

	if (*n == 0)
	{
		return 0;
	}

	/*     Get machine constants */

	eps = dlamch_ ("P", (ftnlen) 1);
	smlnum = dlamch_ ("S", (ftnlen) 1);
	bignum = 1. / smlnum;
	dlabad_ (&smlnum, &bignum);
	smlnum = sqrt (smlnum) / eps;
	bignum = 1. / smlnum;

	/*     Scale A if max element outside range [SMLNUM,BIGNUM] */

	icond = 0;
	anrm = dlange_ ("M", n, n, &a[a_offset], lda, dum, (ftnlen) 1);
	scalea = FALSE_;
	if (anrm > 0. && anrm < smlnum)
	{
		scalea = TRUE_;
		cscale = smlnum;
	}
	else if (anrm > bignum)
	{
		scalea = TRUE_;
		cscale = bignum;
	}
	if (scalea)
	{
		dlascl_ ("G", &c__0, &c__0, &anrm, &cscale, n, n, &a[a_offset], lda, &ierr, (ftnlen) 1);
	}

	/*     Balance the matrix and compute ABNRM */

	dgebal_ (balanc, n, &a[a_offset], lda, ilo, ihi, &scale[1], &ierr, (ftnlen) 1);
	*abnrm = dlange_ ("1", n, n, &a[a_offset], lda, dum, (ftnlen) 1);
	if (scalea)
	{
		dum[0] = *abnrm;
		dlascl_ ("G", &c__0, &c__0, &cscale, &anrm, &c__1, &c__1, dum, &c__1, &ierr, (ftnlen) 1);
		*abnrm = dum[0];
	}

	/*     Reduce to upper Hessenberg form */
	/*     (Workspace: need 2*N, prefer N+N*NB) */

	itau = 1;
	iwrk = itau + *n;
	i__1 = *lwork - iwrk + 1;
	dgehrd_ (n, ilo, ihi, &a[a_offset], lda, &work[itau], &work[iwrk], &i__1, &ierr);

	if (wantvl)
	{

		/*        Want left eigenvectors */
		/*        Copy Householder vectors to VL */

		*(unsigned char *) side = 'L';
		dlacpy_ ("L", n, n, &a[a_offset], lda, &vl[vl_offset], ldvl, (ftnlen) 1);

		/*        Generate orthogonal matrix in VL */
		/*        (Workspace: need 2*N-1, prefer N+(N-1)*NB) */

		i__1 = *lwork - iwrk + 1;
		dorghr_ (n, ilo, ihi, &vl[vl_offset], ldvl, &work[itau], &work[iwrk], &i__1, &ierr);

		/*        Perform QR iteration, accumulating Schur vectors in VL */
		/*        (Workspace: need 1, prefer HSWORK (see comments) ) */

		iwrk = itau;
		i__1 = *lwork - iwrk + 1;
		dhseqr_ ("S", "V", n, ilo, ihi, &a[a_offset], lda, &wr[1], &wi[1], &vl[vl_offset], ldvl, &work[iwrk], &i__1, info, (ftnlen) 1, (ftnlen) 1);

		if (wantvr)
		{

			/*           Want left and right eigenvectors */
			/*           Copy Schur vectors to VR */

			*(unsigned char *) side = 'B';
			dlacpy_ ("F", n, n, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, (ftnlen) 1);
		}

	}
	else if (wantvr)
	{

		/*        Want right eigenvectors */
		/*        Copy Householder vectors to VR */

		*(unsigned char *) side = 'R';
		dlacpy_ ("L", n, n, &a[a_offset], lda, &vr[vr_offset], ldvr, (ftnlen) 1);

		/*        Generate orthogonal matrix in VR */
		/*        (Workspace: need 2*N-1, prefer N+(N-1)*NB) */

		i__1 = *lwork - iwrk + 1;
		dorghr_ (n, ilo, ihi, &vr[vr_offset], ldvr, &work[itau], &work[iwrk], &i__1, &ierr);

		/*        Perform QR iteration, accumulating Schur vectors in VR */
		/*        (Workspace: need 1, prefer HSWORK (see comments) ) */

		iwrk = itau;
		i__1 = *lwork - iwrk + 1;
		dhseqr_ ("S", "V", n, ilo, ihi, &a[a_offset], lda, &wr[1], &wi[1], &vr[vr_offset], ldvr, &work[iwrk], &i__1, info, (ftnlen) 1, (ftnlen) 1);

	}
	else
	{

		/*        Compute eigenvalues only */
		/*        If condition numbers desired, compute Schur form */

		if (wntsnn)
		{
			*(unsigned char *) job = 'E';
		}
		else
		{
			*(unsigned char *) job = 'S';
		}

		/*        (Workspace: need 1, prefer HSWORK (see comments) ) */

		iwrk = itau;
		i__1 = *lwork - iwrk + 1;
		dhseqr_ (job, "N", n, ilo, ihi, &a[a_offset], lda, &wr[1], &wi[1], &vr[vr_offset], ldvr, &work[iwrk], &i__1, info, (ftnlen) 1, (ftnlen) 1);
	}

	/*     If INFO > 0 from DHSEQR, then quit */

	if (*info > 0)
	{
		goto L50;
	}

	if (wantvl || wantvr)
	{

		/*        Compute left and/or right eigenvectors */
		/*        (Workspace: need 3*N) */

		dtrevc_ (side, "B", select, n, &a[a_offset], lda, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, n, &nout, &work[iwrk], &ierr, (ftnlen) 1, (ftnlen) 1);
	}

	/*     Compute condition numbers if desired */
	/*     (Workspace: need N*N+6*N unless SENSE = 'E') */

	if (!wntsnn)
	{
		dtrsna_ (sense, "A", select, n, &a[a_offset], lda, &vl[vl_offset],
					ldvl, &vr[vr_offset], ldvr, &rconde[1], &rcondv[1], n, &nout, &work[iwrk], n, &iwork[1], &icond, (ftnlen) 1, (ftnlen) 1);
	}

	if (wantvl)
	{

		/*        Undo balancing of left eigenvectors */

		dgebak_ (balanc, "L", n, ilo, ihi, &scale[1], n, &vl[vl_offset], ldvl, &ierr, (ftnlen) 1, (ftnlen) 1);

		/*        Normalize left eigenvectors and make largest component real */

		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			if (wi[i__] == 0.)
			{
				scl = 1. / dnrm2_ (n, &vl[i__ * vl_dim1 + 1], &c__1);
				dscal_ (n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
			}
			else if (wi[i__] > 0.)
			{
				d__1 = dnrm2_ (n, &vl[i__ * vl_dim1 + 1], &c__1);
				d__2 = dnrm2_ (n, &vl[(i__ + 1) * vl_dim1 + 1], &c__1);
				scl = 1. / dlapy2_ (&d__1, &d__2);
				dscal_ (n, &scl, &vl[i__ * vl_dim1 + 1], &c__1);
				dscal_ (n, &scl, &vl[(i__ + 1) * vl_dim1 + 1], &c__1);
				i__2 = *n;
				for (k = 1; k <= i__2; ++k)
				{
					/* Computing 2nd power */
					d__1 = vl[k + i__ * vl_dim1];
					/* Computing 2nd power */
					d__2 = vl[k + (i__ + 1) * vl_dim1];
					work[k] = d__1 * d__1 + d__2 * d__2;
					/* L10: */
				}
				k = idamax_ (n, &work[1], &c__1);
				dlartg_ (&vl[k + i__ * vl_dim1], &vl[k + (i__ + 1) * vl_dim1], &cs, &sn, &r__);
				drot_ (n, &vl[i__ * vl_dim1 + 1], &c__1, &vl[(i__ + 1) * vl_dim1 + 1], &c__1, &cs, &sn);
				vl[k + (i__ + 1) * vl_dim1] = 0.;
			}
			/* L20: */
		}
	}

	if (wantvr)
	{

		/*        Undo balancing of right eigenvectors */

		dgebak_ (balanc, "R", n, ilo, ihi, &scale[1], n, &vr[vr_offset], ldvr, &ierr, (ftnlen) 1, (ftnlen) 1);

		/*        Normalize right eigenvectors and make largest component real */

		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			if (wi[i__] == 0.)
			{
				scl = 1. / dnrm2_ (n, &vr[i__ * vr_dim1 + 1], &c__1);
				dscal_ (n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
			}
			else if (wi[i__] > 0.)
			{
				d__1 = dnrm2_ (n, &vr[i__ * vr_dim1 + 1], &c__1);
				d__2 = dnrm2_ (n, &vr[(i__ + 1) * vr_dim1 + 1], &c__1);
				scl = 1. / dlapy2_ (&d__1, &d__2);
				dscal_ (n, &scl, &vr[i__ * vr_dim1 + 1], &c__1);
				dscal_ (n, &scl, &vr[(i__ + 1) * vr_dim1 + 1], &c__1);
				i__2 = *n;
				for (k = 1; k <= i__2; ++k)
				{
					/* Computing 2nd power */
					d__1 = vr[k + i__ * vr_dim1];
					/* Computing 2nd power */
					d__2 = vr[k + (i__ + 1) * vr_dim1];
					work[k] = d__1 * d__1 + d__2 * d__2;
					/* L30: */
				}
				k = idamax_ (n, &work[1], &c__1);
				dlartg_ (&vr[k + i__ * vr_dim1], &vr[k + (i__ + 1) * vr_dim1], &cs, &sn, &r__);
				drot_ (n, &vr[i__ * vr_dim1 + 1], &c__1, &vr[(i__ + 1) * vr_dim1 + 1], &c__1, &cs, &sn);
				vr[k + (i__ + 1) * vr_dim1] = 0.;
			}
			/* L40: */
		}
	}

	/*     Undo scaling if necessary */

 L50:
	if (scalea)
	{
		i__1 = *n - *info;
		/* Computing MAX */
		i__3 = *n - *info;
		i__2 = max (i__3, 1);
		dlascl_ ("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wr[*info + 1], &i__2, &ierr, (ftnlen) 1);
		i__1 = *n - *info;
		/* Computing MAX */
		i__3 = *n - *info;
		i__2 = max (i__3, 1);
		dlascl_ ("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[*info + 1], &i__2, &ierr, (ftnlen) 1);
		if (*info == 0)
		{
			if ((wntsnv || wntsnb) && icond == 0)
			{
				dlascl_ ("G", &c__0, &c__0, &cscale, &anrm, n, &c__1, &rcondv[1], n, &ierr, (ftnlen) 1);
			}
		}
		else
		{
			i__1 = *ilo - 1;
			dlascl_ ("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wr[1], n, &ierr, (ftnlen) 1);
			i__1 = *ilo - 1;
			dlascl_ ("G", &c__0, &c__0, &cscale, &anrm, &i__1, &c__1, &wi[1], n, &ierr, (ftnlen) 1);
		}
	}

	work[1] = (doublereal) maxwrk;
	return 0;

	/*     End of DGEEVX */

}	/* dgeevx_ */

static int
dgehd2_ (integer * n, integer * ilo, integer * ihi, doublereal * a, integer * lda, doublereal * tau, doublereal * work, integer * info)
{
	/* System generated locals */
	integer a_dim1, a_offset, i__1, i__2, i__3;

	/* Local variables */
	static integer i__;
	static doublereal aii;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     October 31, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DGEHD2 reduces a real general matrix A to upper Hessenberg form H by */
	/*  an orthogonal similarity transformation:  Q' * A * Q = H . */

	/*  Arguments */
	/*  ========= */

	/*  N       (input) INTEGER */
	/*          The order of the matrix A.  N >= 0. */

	/*  ILO     (input) INTEGER */
	/*  IHI     (input) INTEGER */
	/*          It is assumed that A is already upper triangular in rows */
	/*          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally */
	/*          set by a previous call to DGEBAL; otherwise they should be */
	/*          set to 1 and N respectively. See Further Details. */
	/*          1 <= ILO <= IHI <= max(1,N). */

	/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          On entry, the n by n general matrix to be reduced. */
	/*          On exit, the upper triangle and the first subdiagonal of A */
	/*          are overwritten with the upper Hessenberg matrix H, and the */
	/*          elements below the first subdiagonal, with the array TAU, */
	/*          represent the orthogonal matrix Q as a product of elementary */
	/*          reflectors. See Further Details. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(1,N). */

	/*  TAU     (output) DOUBLE PRECISION array, dimension (N-1) */
	/*          The scalar factors of the elementary reflectors (see Further */
	/*          Details). */

	/*  WORK    (workspace) DOUBLE PRECISION array, dimension (N) */

	/*  INFO    (output) INTEGER */
	/*          = 0:  successful exit. */
	/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */

	/*  Further Details */
	/*  =============== */

	/*  The matrix Q is represented as a product of (ihi-ilo) elementary */
	/*  reflectors */

	/*     Q = H(ilo) H(ilo+1) . . . H(ihi-1). */

	/*  Each H(i) has the form */

	/*     H(i) = I - tau * v * v' */

	/*  where tau is a real scalar, and v is a real vector with */
	/*  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on */
	/*  exit in A(i+2:ihi,i), and tau in TAU(i). */

	/*  The contents of A are illustrated by the following example, with */
	/*  n = 7, ilo = 2 and ihi = 6: */

	/*  on entry,                        on exit, */

	/*  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a ) */
	/*  (     a   a   a   a   a   a )    (      a   h   h   h   h   a ) */
	/*  (     a   a   a   a   a   a )    (      h   h   h   h   h   h ) */
	/*  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h ) */
	/*  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h ) */
	/*  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h ) */
	/*  (                         a )    (                          a ) */

	/*  where a denotes an element of the original matrix A, h denotes a */
	/*  modified element of the upper Hessenberg matrix H, and vi denotes an */
	/*  element of the vector defining H(i). */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Test the input parameters */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--tau;
	--work;

	/* Function Body */
	*info = 0;
	if (*n < 0)
	{
		*info = -1;
	}
	else if (*ilo < 1 || *ilo > max (1, *n))
	{
		*info = -2;
	}
	else if (*ihi < min (*ilo, *n) || *ihi > *n)
	{
		*info = -3;
	}
	else if (*lda < max (1, *n))
	{
		*info = -5;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DGEHD2", &i__1, (ftnlen) 6);
		return 0;
	}

	i__1 = *ihi - 1;
	for (i__ = *ilo; i__ <= i__1; ++i__)
	{

		/*        Compute elementary reflector H(i) to annihilate A(i+2:ihi,i) */

		i__2 = *ihi - i__;
		/* Computing MIN */
		i__3 = i__ + 2;
		dlarfg_ (&i__2, &a[i__ + 1 + i__ * a_dim1], &a[min (i__3, *n) + i__ * a_dim1], &c__1, &tau[i__]);
		aii = a[i__ + 1 + i__ * a_dim1];
		a[i__ + 1 + i__ * a_dim1] = 1.;

		/*        Apply H(i) to A(1:ihi,i+1:ihi) from the right */

		i__2 = *ihi - i__;
		dlarf_ ("Right", ihi, &i__2, &a[i__ + 1 + i__ * a_dim1], &c__1, &tau[i__], &a[(i__ + 1) * a_dim1 + 1], lda, &work[1], (ftnlen) 5);

		/*        Apply H(i) to A(i+1:ihi,i+1:n) from the left */

		i__2 = *ihi - i__;
		i__3 = *n - i__;
		dlarf_ ("Left", &i__2, &i__3, &a[i__ + 1 + i__ * a_dim1], &c__1, &tau[i__], &a[i__ + 1 + (i__ + 1) * a_dim1], lda, &work[1], (ftnlen) 4);

		a[i__ + 1 + i__ * a_dim1] = aii;
		/* L10: */
	}

	return 0;

	/*     End of DGEHD2 */

}	/* dgehd2_ */

static int
dgehrd_ (integer * n, integer * ilo, integer * ihi, doublereal * a, integer * lda, doublereal * tau, doublereal * work, integer * lwork, integer * info)
{
	/* System generated locals */
	integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

	/* Local variables */
	static integer i__;
	static doublereal t[4160] /* was [65][64] */ ;
	static integer ib;
	static doublereal ei;
	static integer nb, nh, nx, iws;
	static integer nbmin, iinfo;
	static integer ldwork, lwkopt;
	static logical lquery;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     June 30, 1999 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DGEHRD reduces a real general matrix A to upper Hessenberg form H by */
	/*  an orthogonal similarity transformation:  Q' * A * Q = H . */

	/*  Arguments */
	/*  ========= */

	/*  N       (input) INTEGER */
	/*          The order of the matrix A.  N >= 0. */

	/*  ILO     (input) INTEGER */
	/*  IHI     (input) INTEGER */
	/*          It is assumed that A is already upper triangular in rows */
	/*          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally */
	/*          set by a previous call to DGEBAL; otherwise they should be */
	/*          set to 1 and N respectively. See Further Details. */
	/*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */

	/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          On entry, the N-by-N general matrix to be reduced. */
	/*          On exit, the upper triangle and the first subdiagonal of A */
	/*          are overwritten with the upper Hessenberg matrix H, and the */
	/*          elements below the first subdiagonal, with the array TAU, */
	/*          represent the orthogonal matrix Q as a product of elementary */
	/*          reflectors. See Further Details. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(1,N). */

	/*  TAU     (output) DOUBLE PRECISION array, dimension (N-1) */
	/*          The scalar factors of the elementary reflectors (see Further */
	/*          Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to */
	/*          zero. */

	/*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) */
	/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

	/*  LWORK   (input) INTEGER */
	/*          The length of the array WORK.  LWORK >= max(1,N). */
	/*          For optimum performance LWORK >= N*NB, where NB is the */
	/*          optimal blocksize. */

	/*          If LWORK = -1, then a workspace query is assumed; the routine */
	/*          only calculates the optimal size of the WORK array, returns */
	/*          this value as the first entry of the WORK array, and no error */
	/*          message related to LWORK is issued by XERBLA. */

	/*  INFO    (output) INTEGER */
	/*          = 0:  successful exit */
	/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */

	/*  Further Details */
	/*  =============== */

	/*  The matrix Q is represented as a product of (ihi-ilo) elementary */
	/*  reflectors */

	/*     Q = H(ilo) H(ilo+1) . . . H(ihi-1). */

	/*  Each H(i) has the form */

	/*     H(i) = I - tau * v * v' */

	/*  where tau is a real scalar, and v is a real vector with */
	/*  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on */
	/*  exit in A(i+2:ihi,i), and tau in TAU(i). */

	/*  The contents of A are illustrated by the following example, with */
	/*  n = 7, ilo = 2 and ihi = 6: */

	/*  on entry,                        on exit, */

	/*  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a ) */
	/*  (     a   a   a   a   a   a )    (      a   h   h   h   h   a ) */
	/*  (     a   a   a   a   a   a )    (      h   h   h   h   h   h ) */
	/*  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h ) */
	/*  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h ) */
	/*  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h ) */
	/*  (                         a )    (                          a ) */

	/*  where a denotes an element of the original matrix A, h denotes a */
	/*  modified element of the upper Hessenberg matrix H, and vi denotes an */
	/*  element of the vector defining H(i). */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Local Arrays .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Test the input parameters */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--tau;
	--work;

	/* Function Body */
	*info = 0;
	/* Computing MIN */
	i__1 = 64, i__2 = ilaenv_ (&c__1, "DGEHRD", " ", n, ilo, ihi, &c_n1, (ftnlen) 6, (ftnlen) 1);
	nb = min (i__1, i__2);
	lwkopt = *n * nb;
	work[1] = (doublereal) lwkopt;
	lquery = *lwork == -1;
	if (*n < 0)
	{
		*info = -1;
	}
	else if (*ilo < 1 || *ilo > max (1, *n))
	{
		*info = -2;
	}
	else if (*ihi < min (*ilo, *n) || *ihi > *n)
	{
		*info = -3;
	}
	else if (*lda < max (1, *n))
	{
		*info = -5;
	}
	else if (*lwork < max (1, *n) && !lquery)
	{
		*info = -8;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DGEHRD", &i__1, (ftnlen) 6);
		return 0;
	}
	else if (lquery)
	{
		return 0;
	}

	/*     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero */

	i__1 = *ilo - 1;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		tau[i__] = 0.;
		/* L10: */
	}
	i__1 = *n - 1;
	for (i__ = max (1, *ihi); i__ <= i__1; ++i__)
	{
		tau[i__] = 0.;
		/* L20: */
	}

	/*     Quick return if possible */

	nh = *ihi - *ilo + 1;
	if (nh <= 1)
	{
		work[1] = 1.;
		return 0;
	}

	/*     Determine the block size. */

	/* Computing MIN */
	i__1 = 64, i__2 = ilaenv_ (&c__1, "DGEHRD", " ", n, ilo, ihi, &c_n1, (ftnlen) 6, (ftnlen) 1);
	nb = min (i__1, i__2);
	nbmin = 2;
	iws = 1;
	if (nb > 1 && nb < nh)
	{

		/*        Determine when to cross over from blocked to unblocked code */
		/*        (last block is always handled by unblocked code). */

		/* Computing MAX */
		i__1 = nb, i__2 = ilaenv_ (&c__3, "DGEHRD", " ", n, ilo, ihi, &c_n1, (ftnlen) 6, (ftnlen) 1);
		nx = max (i__1, i__2);
		if (nx < nh)
		{

			/*           Determine if workspace is large enough for blocked code. */

			iws = *n * nb;
			if (*lwork < iws)
			{

				/*              Not enough workspace to use optimal NB:  determine the */
				/*              minimum value of NB, and reduce NB or force use of */
				/*              unblocked code. */

				/* Computing MAX */
				i__1 = 2, i__2 = ilaenv_ (&c__2, "DGEHRD", " ", n, ilo, ihi, &c_n1, (ftnlen) 6, (ftnlen) 1);
				nbmin = max (i__1, i__2);
				if (*lwork >= *n * nbmin)
				{
					nb = *lwork / *n;
				}
				else
				{
					nb = 1;
				}
			}
		}
	}
	ldwork = *n;

	if (nb < nbmin || nb >= nh)
	{

		/*        Use unblocked code below */

		i__ = *ilo;

	}
	else
	{

		/*        Use blocked code */

		i__1 = *ihi - 1 - nx;
		i__2 = nb;
		for (i__ = *ilo; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
		{
			/* Computing MIN */
			i__3 = nb, i__4 = *ihi - i__;
			ib = min (i__3, i__4);

			/*           Reduce columns i:i+ib-1 to Hessenberg form, returning the */
			/*           matrices V and T of the block reflector H = I - V*T*V' */
			/*           which performs the reduction, and also the matrix Y = A*V*T */

			dlahrd_ (ihi, &i__, &ib, &a[i__ * a_dim1 + 1], lda, &tau[i__], t, &c__65, &work[1], &ldwork);

			/*           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the */
			/*           right, computing  A := A - Y * V'. V(i+ib,ib-1) must be set */
			/*           to 1. */

			ei = a[i__ + ib + (i__ + ib - 1) * a_dim1];
			a[i__ + ib + (i__ + ib - 1) * a_dim1] = 1.;
			i__3 = *ihi - i__ - ib + 1;
			dgemm_ ("No transpose", "Transpose", ihi, &i__3, &ib, &c_b347, &work[1], &ldwork, &a[i__ + ib + i__ * a_dim1], lda, &c_b348,
					  &a[(i__ + ib) * a_dim1 + 1], lda, (ftnlen) 12, (ftnlen) 9);
			a[i__ + ib + (i__ + ib - 1) * a_dim1] = ei;

			/*           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the */
			/*           left */

			i__3 = *ihi - i__;
			i__4 = *n - i__ - ib + 1;
			dlarfb_ ("Left", "Transpose", "Forward", "Columnwise", &i__3, &i__4, &ib, &a[i__ + 1 + i__ * a_dim1], lda, t, &c__65,
						&a[i__ + 1 + (i__ + ib) * a_dim1], lda, &work[1], &ldwork, (ftnlen) 4, (ftnlen) 9, (ftnlen) 7, (ftnlen) 10);
			/* L30: */
		}
	}

	/*     Use unblocked code to reduce the rest of the matrix */

	dgehd2_ (n, &i__, ihi, &a[a_offset], lda, &tau[1], &work[1], &iinfo);
	work[1] = (doublereal) iws;

	return 0;

	/*     End of DGEHRD */

}	/* dgehrd_ */

static int
dgeqr2_ (integer * m, integer * n, doublereal * a, integer * lda, doublereal * tau, doublereal * work, integer * info)
{
	/* System generated locals */
	integer a_dim1, a_offset, i__1, i__2, i__3;

	/* Local variables */
	static integer i__, k;
	static doublereal aii;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     February 29, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DGEQR2 computes a QR factorization of a real m by n matrix A: */
	/*  A = Q * R. */

	/*  Arguments */
	/*  ========= */

	/*  M       (input) INTEGER */
	/*          The number of rows of the matrix A.  M >= 0. */

	/*  N       (input) INTEGER */
	/*          The number of columns of the matrix A.  N >= 0. */

	/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          On entry, the m by n matrix A. */
	/*          On exit, the elements on and above the diagonal of the array */
	/*          contain the min(m,n) by n upper trapezoidal matrix R (R is */
	/*          upper triangular if m >= n); the elements below the diagonal, */
	/*          with the array TAU, represent the orthogonal matrix Q as a */
	/*          product of elementary reflectors (see Further Details). */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(1,M). */

	/*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N)) */
	/*          The scalar factors of the elementary reflectors (see Further */
	/*          Details). */

	/*  WORK    (workspace) DOUBLE PRECISION array, dimension (N) */

	/*  INFO    (output) INTEGER */
	/*          = 0: successful exit */
	/*          < 0: if INFO = -i, the i-th argument had an illegal value */

	/*  Further Details */
	/*  =============== */

	/*  The matrix Q is represented as a product of elementary reflectors */

	/*     Q = H(1) H(2) . . . H(k), where k = min(m,n). */

	/*  Each H(i) has the form */

	/*     H(i) = I - tau * v * v' */

	/*  where tau is a real scalar, and v is a real vector with */
	/*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i), */
	/*  and tau in TAU(i). */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Test the input arguments */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--tau;
	--work;

	/* Function Body */
	*info = 0;
	if (*m < 0)
	{
		*info = -1;
	}
	else if (*n < 0)
	{
		*info = -2;
	}
	else if (*lda < max (1, *m))
	{
		*info = -4;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DGEQR2", &i__1, (ftnlen) 6);
		return 0;
	}

	k = min (*m, *n);

	i__1 = k;
	for (i__ = 1; i__ <= i__1; ++i__)
	{

		/*        Generate elementary reflector H(i) to annihilate A(i+1:m,i) */

		i__2 = *m - i__ + 1;
		/* Computing MIN */
		i__3 = i__ + 1;
		dlarfg_ (&i__2, &a[i__ + i__ * a_dim1], &a[min (i__3, *m) + i__ * a_dim1], &c__1, &tau[i__]);
		if (i__ < *n)
		{

			/*           Apply H(i) to A(i:m,i+1:n) from the left */

			aii = a[i__ + i__ * a_dim1];
			a[i__ + i__ * a_dim1] = 1.;
			i__2 = *m - i__ + 1;
			i__3 = *n - i__;
			dlarf_ ("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &tau[i__], &a[i__ + (i__ + 1) * a_dim1], lda, &work[1], (ftnlen) 4);
			a[i__ + i__ * a_dim1] = aii;
		}
		/* L10: */
	}
	return 0;

	/*     End of DGEQR2 */

}	/* dgeqr2_ */

static int
dgerfs_ (char *trans, integer * n, integer * nrhs,
			doublereal * a, integer * lda, doublereal * af, integer * ldaf, integer *
			ipiv, doublereal * b, integer * ldb, doublereal * x, integer * ldx,
			doublereal * ferr, doublereal * berr, doublereal * work, integer * iwork, integer * info, ftnlen trans_len)
{
	/* System generated locals */
	integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, x_dim1, x_offset, i__1, i__2, i__3;
	doublereal d__1, d__2, d__3;

	/* Local variables */
	static integer i__, j, k;
	static doublereal s, xk;
	static integer nz;
	static doublereal eps;
	static integer kase;
	static doublereal safe1, safe2;
	static integer count;
	static doublereal safmin;
	static logical notran;
	static char transt[1];
	static doublereal lstres;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     September 30, 1994 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DGERFS improves the computed solution to a system of linear */
	/*  equations and provides error bounds and backward error estimates for */
	/*  the solution. */

	/*  Arguments */
	/*  ========= */

	/*  TRANS   (input) CHARACTER*1 */
	/*          Specifies the form of the system of equations: */
	/*          = 'N':  A * X = B     (No transpose) */
	/*          = 'T':  A**T * X = B  (Transpose) */
	/*          = 'C':  A**H * X = B  (Conjugate transpose = Transpose) */

	/*  N       (input) INTEGER */
	/*          The order of the matrix A.  N >= 0. */

	/*  NRHS    (input) INTEGER */
	/*          The number of right hand sides, i.e., the number of columns */
	/*          of the matrices B and X.  NRHS >= 0. */

	/*  A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          The original N-by-N matrix A. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(1,N). */

	/*  AF      (input) DOUBLE PRECISION array, dimension (LDAF,N) */
	/*          The factors L and U from the factorization A = P*L*U */
	/*          as computed by DGETRF. */

	/*  LDAF    (input) INTEGER */
	/*          The leading dimension of the array AF.  LDAF >= max(1,N). */

	/*  IPIV    (input) INTEGER array, dimension (N) */
	/*          The pivot indices from DGETRF; for 1<=i<=N, row i of the */
	/*          matrix was interchanged with row IPIV(i). */

	/*  B       (input) DOUBLE PRECISION array, dimension (LDB,NRHS) */
	/*          The right hand side matrix B. */

	/*  LDB     (input) INTEGER */
	/*          The leading dimension of the array B.  LDB >= max(1,N). */

	/*  X       (input/output) DOUBLE PRECISION array, dimension (LDX,NRHS) */
	/*          On entry, the solution matrix X, as computed by DGETRS. */
	/*          On exit, the improved solution matrix X. */

	/*  LDX     (input) INTEGER */
	/*          The leading dimension of the array X.  LDX >= max(1,N). */

	/*  FERR    (output) DOUBLE PRECISION array, dimension (NRHS) */
	/*          The estimated forward error bound for each solution vector */
	/*          X(j) (the j-th column of the solution matrix X). */
	/*          If XTRUE is the true solution corresponding to X(j), FERR(j) */
	/*          is an estimated upper bound for the magnitude of the largest */
	/*          element in (X(j) - XTRUE) divided by the magnitude of the */
	/*          largest element in X(j).  The estimate is as reliable as */
	/*          the estimate for RCOND, and is almost always a slight */
	/*          overestimate of the true error. */

	/*  BERR    (output) DOUBLE PRECISION array, dimension (NRHS) */
	/*          The componentwise relative backward error of each solution */
	/*          vector X(j) (i.e., the smallest relative change in */
	/*          any element of A or B that makes X(j) an exact solution). */

	/*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N) */

	/*  IWORK   (workspace) INTEGER array, dimension (N) */

	/*  INFO    (output) INTEGER */
	/*          = 0:  successful exit */
	/*          < 0:  if INFO = -i, the i-th argument had an illegal value */

	/*  Internal Parameters */
	/*  =================== */

	/*  ITMAX is the maximum number of steps of iterative refinement. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Test the input parameters. */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	af_dim1 = *ldaf;
	af_offset = 1 + af_dim1;
	af -= af_offset;
	--ipiv;
	b_dim1 = *ldb;
	b_offset = 1 + b_dim1;
	b -= b_offset;
	x_dim1 = *ldx;
	x_offset = 1 + x_dim1;
	x -= x_offset;
	--ferr;
	--berr;
	--work;
	--iwork;

	/* Function Body */
	*info = 0;
	notran = lsame_ (trans, "N", (ftnlen) 1, (ftnlen) 1);
	if (!notran && !lsame_ (trans, "T", (ftnlen) 1, (ftnlen) 1) && !lsame_ (trans, "C", (ftnlen) 1, (ftnlen) 1))
	{
		*info = -1;
	}
	else if (*n < 0)
	{
		*info = -2;
	}
	else if (*nrhs < 0)
	{
		*info = -3;
	}
	else if (*lda < max (1, *n))
	{
		*info = -5;
	}
	else if (*ldaf < max (1, *n))
	{
		*info = -7;
	}
	else if (*ldb < max (1, *n))
	{
		*info = -10;
	}
	else if (*ldx < max (1, *n))
	{
		*info = -12;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DGERFS", &i__1, (ftnlen) 6);
		return 0;
	}

	/*     Quick return if possible */

	if (*n == 0 || *nrhs == 0)
	{
		i__1 = *nrhs;
		for (j = 1; j <= i__1; ++j)
		{
			ferr[j] = 0.;
			berr[j] = 0.;
			/* L10: */
		}
		return 0;
	}

	if (notran)
	{
		*(unsigned char *) transt = 'T';
	}
	else
	{
		*(unsigned char *) transt = 'N';
	}

	/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

	nz = *n + 1;
	eps = dlamch_ ("Epsilon", (ftnlen) 7);
	safmin = dlamch_ ("Safe minimum", (ftnlen) 12);
	safe1 = nz * safmin;
	safe2 = safe1 / eps;

	/*     Do for each right hand side */

	i__1 = *nrhs;
	for (j = 1; j <= i__1; ++j)
	{

		count = 1;
		lstres = 3.;
	 L20:

		/*        Loop until stopping criterion is satisfied. */

		/*        Compute residual R = B - op(A) * X, */
		/*        where op(A) = A, A**T, or A**H, depending on TRANS. */

		dcopy_ (n, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);
		dgemv_ (trans, n, n, &c_b347, &a[a_offset], lda, &x[j * x_dim1 + 1], &c__1, &c_b348, &work[*n + 1], &c__1, (ftnlen) 1);

		/*        Compute componentwise relative backward error from formula */

		/*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) ) */

		/*        where abs(Z) is the componentwise absolute value of the matrix */
		/*        or vector Z.  If the i-th component of the denominator is less */
		/*        than SAFE2, then SAFE1 is added to the i-th components of the */
		/*        numerator and denominator before dividing. */

		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__)
		{
			work[i__] = (d__1 = b[i__ + j * b_dim1], abs (d__1));
			/* L30: */
		}

		/*        Compute abs(op(A))*abs(X) + abs(B). */

		if (notran)
		{
			i__2 = *n;
			for (k = 1; k <= i__2; ++k)
			{
				xk = (d__1 = x[k + j * x_dim1], abs (d__1));
				i__3 = *n;
				for (i__ = 1; i__ <= i__3; ++i__)
				{
					work[i__] += (d__1 = a[i__ + k * a_dim1], abs (d__1)) * xk;
					/* L40: */
				}
				/* L50: */
			}
		}
		else
		{
			i__2 = *n;
			for (k = 1; k <= i__2; ++k)
			{
				s = 0.;
				i__3 = *n;
				for (i__ = 1; i__ <= i__3; ++i__)
				{
					s += (d__1 = a[i__ + k * a_dim1], abs (d__1)) * (d__2 = x[i__ + j * x_dim1], abs (d__2));
					/* L60: */
				}
				work[k] += s;
				/* L70: */
			}
		}
		s = 0.;
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__)
		{
			if (work[i__] > safe2)
			{
				/* Computing MAX */
				d__2 = s, d__3 = (d__1 = work[*n + i__], abs (d__1)) / work[i__];
				s = max (d__2, d__3);
			}
			else
			{
				/* Computing MAX */
				d__2 = s, d__3 = ((d__1 = work[*n + i__], abs (d__1)) + safe1) / (work[i__] + safe1);
				s = max (d__2, d__3);
			}
			/* L80: */
		}
		berr[j] = s;

		/*        Test stopping criterion. Continue iterating if */
		/*           1) The residual BERR(J) is larger than machine epsilon, and */
		/*           2) BERR(J) decreased by at least a factor of 2 during the */
		/*              last iteration, and */
		/*           3) At most ITMAX iterations tried. */

		if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5)
		{

			/*           Update solution and try again. */

			dgetrs_ (trans, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[*n + 1], n, info, (ftnlen) 1);
			daxpy_ (n, &c_b348, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1);
			lstres = berr[j];
			++count;
			goto L20;
		}

		/*        Bound error from formula */

		/*        norm(X - XTRUE) / norm(X) .le. FERR = */
		/*        norm( abs(inv(op(A)))* */
		/*           ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X) */

		/*        where */
		/*          norm(Z) is the magnitude of the largest component of Z */
		/*          inv(op(A)) is the inverse of op(A) */
		/*          abs(Z) is the componentwise absolute value of the matrix or */
		/*             vector Z */
		/*          NZ is the maximum number of nonzeros in any row of A, plus 1 */
		/*          EPS is machine epsilon */

		/*        The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B)) */
		/*        is incremented by SAFE1 if the i-th component of */
		/*        abs(op(A))*abs(X) + abs(B) is less than SAFE2. */

		/*        Use DLACON to estimate the infinity-norm of the matrix */
		/*           inv(op(A)) * diag(W), */
		/*        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) */

		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__)
		{
			if (work[i__] > safe2)
			{
				work[i__] = (d__1 = work[*n + i__], abs (d__1)) + nz * eps * work[i__];
			}
			else
			{
				work[i__] = (d__1 = work[*n + i__], abs (d__1)) + nz * eps * work[i__] + safe1;
			}
			/* L90: */
		}

		kase = 0;
	 L100:
		dlacon_ (n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &kase);
		if (kase != 0)
		{
			if (kase == 1)
			{

				/*              Multiply by diag(W)*inv(op(A)**T). */

				dgetrs_ (transt, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[*n + 1], n, info, (ftnlen) 1);
				i__2 = *n;
				for (i__ = 1; i__ <= i__2; ++i__)
				{
					work[*n + i__] = work[i__] * work[*n + i__];
					/* L110: */
				}
			}
			else
			{

				/*              Multiply by inv(op(A))*diag(W). */

				i__2 = *n;
				for (i__ = 1; i__ <= i__2; ++i__)
				{
					work[*n + i__] = work[i__] * work[*n + i__];
					/* L120: */
				}
				dgetrs_ (trans, n, &c__1, &af[af_offset], ldaf, &ipiv[1], &work[*n + 1], n, info, (ftnlen) 1);
			}
			goto L100;
		}

		/*        Normalize error. */

		lstres = 0.;
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__)
		{
			/* Computing MAX */
			d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs (d__1));
			lstres = max (d__2, d__3);
			/* L130: */
		}
		if (lstres != 0.)
		{
			ferr[j] /= lstres;
		}

		/* L140: */
	}

	return 0;

	/*     End of DGERFS */

}	/* dgerfs_ */

int
dgesv_ (integer * n, integer * nrhs, doublereal * a, integer * lda, integer * ipiv, doublereal * b, integer * ldb, integer * info)
{
	/* System generated locals */
	integer a_dim1, a_offset, b_dim1, b_offset, i__1;

	/* Local variables */


	/*  -- LAPACK driver routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     March 31, 1993 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DGESV computes the solution to a real system of linear equations */
	/*     A * X = B, */
	/*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices. */

	/*  The LU decomposition with partial pivoting and row interchanges is */
	/*  used to factor A as */
	/*     A = P * L * U, */
	/*  where P is a permutation matrix, L is unit lower triangular, and U is */
	/*  upper triangular.  The factored form of A is then used to solve the */
	/*  system of equations A * X = B. */

	/*  Arguments */
	/*  ========= */

	/*  N       (input) INTEGER */
	/*          The number of linear equations, i.e., the order of the */
	/*          matrix A.  N >= 0. */

	/*  NRHS    (input) INTEGER */
	/*          The number of right hand sides, i.e., the number of columns */
	/*          of the matrix B.  NRHS >= 0. */

	/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          On entry, the N-by-N coefficient matrix A. */
	/*          On exit, the factors L and U from the factorization */
	/*          A = P*L*U; the unit diagonal elements of L are not stored. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(1,N). */

	/*  IPIV    (output) INTEGER array, dimension (N) */
	/*          The pivot indices that define the permutation matrix P; */
	/*          row i of the matrix was interchanged with row IPIV(i). */

	/*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS) */
	/*          On entry, the N-by-NRHS matrix of right hand side matrix B. */
	/*          On exit, if INFO = 0, the N-by-NRHS solution matrix X. */

	/*  LDB     (input) INTEGER */
	/*          The leading dimension of the array B.  LDB >= max(1,N). */

	/*  INFO    (output) INTEGER */
	/*          = 0:  successful exit */
	/*          < 0:  if INFO = -i, the i-th argument had an illegal value */
	/*          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization */
	/*                has been completed, but the factor U is exactly */
	/*                singular, so the solution could not be computed. */

	/*  ===================================================================== */

	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Test the input parameters. */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--ipiv;
	b_dim1 = *ldb;
	b_offset = 1 + b_dim1;
	b -= b_offset;

	/* Function Body */
	*info = 0;
	if (*n < 0)
	{
		*info = -1;
	}
	else if (*nrhs < 0)
	{
		*info = -2;
	}
	else if (*lda < max (1, *n))
	{
		*info = -4;
	}
	else if (*ldb < max (1, *n))
	{
		*info = -7;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DGESV ", &i__1, (ftnlen) 6);
		return 0;
	}

	/*     Compute the LU factorization of A. */

	dgetrf_ (n, n, &a[a_offset], lda, &ipiv[1], info);
	if (*info == 0)
	{

		/*        Solve the system A*X = B, overwriting B with X. */

		dgetrs_ ("No transpose", n, nrhs, &a[a_offset], lda, &ipiv[1], &b[b_offset], ldb, info, (ftnlen) 12);
	}
	return 0;

	/*     End of DGESV */

}	/* dgesv_ */

int
dgesvx_ (char *fact, char *trans, integer * n, integer *
			nrhs, doublereal * a, integer * lda, doublereal * af, integer * ldaf,
			integer * ipiv, char *equed, doublereal * r__, doublereal * c__,
			doublereal * b, integer * ldb, doublereal * x, integer * ldx, doublereal *
			rcond, doublereal * ferr, doublereal * berr, doublereal * work, integer * iwork, integer * info, ftnlen fact_len, ftnlen trans_len, ftnlen equed_len)
{
	/* System generated locals */
	integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, x_dim1, x_offset, i__1, i__2;
	doublereal d__1, d__2;

	/* Local variables */
	static integer i__, j;
	static doublereal amax;
	static char norm[1];
	static doublereal rcmin, rcmax, anorm;
	static logical equil;
	static doublereal colcnd;
	static logical nofact;
	static doublereal bignum;
	static integer infequ;
	static logical colequ;
	static doublereal rowcnd;
	static logical notran;
	static doublereal smlnum;
	static logical rowequ;
	static doublereal rpvgrw;


	/*  -- LAPACK driver routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     June 30, 1999 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DGESVX uses the LU factorization to compute the solution to a real */
	/*  system of linear equations */
	/*     A * X = B, */
	/*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices. */

	/*  Error bounds on the solution and a condition estimate are also */
	/*  provided. */

	/*  Description */
	/*  =========== */

	/*  The following steps are performed: */

	/*  1. If FACT = 'E', real scaling factors are computed to equilibrate */
	/*     the system: */
	/*        TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B */
	/*        TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B */
	/*        TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B */
	/*     Whether or not the system will be equilibrated depends on the */
	/*     scaling of the matrix A, but if equilibration is used, A is */
	/*     overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N') */
	/*     or diag(C)*B (if TRANS = 'T' or 'C'). */

	/*  2. If FACT = 'N' or 'E', the LU decomposition is used to factor the */
	/*     matrix A (after equilibration if FACT = 'E') as */
	/*        A = P * L * U, */
	/*     where P is a permutation matrix, L is a unit lower triangular */
	/*     matrix, and U is upper triangular. */

	/*  3. If some U(i,i)=0, so that U is exactly singular, then the routine */
	/*     returns with INFO = i. Otherwise, the factored form of A is used */
	/*     to estimate the condition number of the matrix A.  If the */
	/*     reciprocal of the condition number is less than machine precision, */
	/*     INFO = N+1 is returned as a warning, but the routine still goes on */
	/*     to solve for X and compute error bounds as described below. */

	/*  4. The system of equations is solved for X using the factored form */
	/*     of A. */

	/*  5. Iterative refinement is applied to improve the computed solution */
	/*     matrix and calculate error bounds and backward error estimates */
	/*     for it. */

	/*  6. If equilibration was used, the matrix X is premultiplied by */
	/*     diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so */
	/*     that it solves the original system before equilibration. */

	/*  Arguments */
	/*  ========= */

	/*  FACT    (input) CHARACTER*1 */
	/*          Specifies whether or not the factored form of the matrix A is */
	/*          supplied on entry, and if not, whether the matrix A should be */
	/*          equilibrated before it is factored. */
	/*          = 'F':  On entry, AF and IPIV contain the factored form of A. */
	/*                  If EQUED is not 'N', the matrix A has been */
	/*                  equilibrated with scaling factors given by R and C. */
	/*                  A, AF, and IPIV are not modified. */
	/*          = 'N':  The matrix A will be copied to AF and factored. */
	/*          = 'E':  The matrix A will be equilibrated if necessary, then */
	/*                  copied to AF and factored. */

	/*  TRANS   (input) CHARACTER*1 */
	/*          Specifies the form of the system of equations: */
	/*          = 'N':  A * X = B     (No transpose) */
	/*          = 'T':  A**T * X = B  (Transpose) */
	/*          = 'C':  A**H * X = B  (Transpose) */

	/*  N       (input) INTEGER */
	/*          The number of linear equations, i.e., the order of the */
	/*          matrix A.  N >= 0. */

	/*  NRHS    (input) INTEGER */
	/*          The number of right hand sides, i.e., the number of columns */
	/*          of the matrices B and X.  NRHS >= 0. */

	/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          On entry, the N-by-N matrix A.  If FACT = 'F' and EQUED is */
	/*          not 'N', then A must have been equilibrated by the scaling */
	/*          factors in R and/or C.  A is not modified if FACT = 'F' or */
	/*          'N', or if FACT = 'E' and EQUED = 'N' on exit. */

	/*          On exit, if EQUED .ne. 'N', A is scaled as follows: */
	/*          EQUED = 'R':  A := diag(R) * A */
	/*          EQUED = 'C':  A := A * diag(C) */
	/*          EQUED = 'B':  A := diag(R) * A * diag(C). */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(1,N). */

	/*  AF      (input or output) DOUBLE PRECISION array, dimension (LDAF,N) */
	/*          If FACT = 'F', then AF is an input argument and on entry */
	/*          contains the factors L and U from the factorization */
	/*          A = P*L*U as computed by DGETRF.  If EQUED .ne. 'N', then */
	/*          AF is the factored form of the equilibrated matrix A. */

	/*          If FACT = 'N', then AF is an output argument and on exit */
	/*          returns the factors L and U from the factorization A = P*L*U */
	/*          of the original matrix A. */

	/*          If FACT = 'E', then AF is an output argument and on exit */
	/*          returns the factors L and U from the factorization A = P*L*U */
	/*          of the equilibrated matrix A (see the description of A for */
	/*          the form of the equilibrated matrix). */

	/*  LDAF    (input) INTEGER */
	/*          The leading dimension of the array AF.  LDAF >= max(1,N). */

	/*  IPIV    (input or output) INTEGER array, dimension (N) */
	/*          If FACT = 'F', then IPIV is an input argument and on entry */
	/*          contains the pivot indices from the factorization A = P*L*U */
	/*          as computed by DGETRF; row i of the matrix was interchanged */
	/*          with row IPIV(i). */

	/*          If FACT = 'N', then IPIV is an output argument and on exit */
	/*          contains the pivot indices from the factorization A = P*L*U */
	/*          of the original matrix A. */

	/*          If FACT = 'E', then IPIV is an output argument and on exit */
	/*          contains the pivot indices from the factorization A = P*L*U */
	/*          of the equilibrated matrix A. */

	/*  EQUED   (input or output) CHARACTER*1 */
	/*          Specifies the form of equilibration that was done. */
	/*          = 'N':  No equilibration (always true if FACT = 'N'). */
	/*          = 'R':  Row equilibration, i.e., A has been premultiplied by */
	/*                  diag(R). */
	/*          = 'C':  Column equilibration, i.e., A has been postmultiplied */
	/*                  by diag(C). */
	/*          = 'B':  Both row and column equilibration, i.e., A has been */
	/*                  replaced by diag(R) * A * diag(C). */
	/*          EQUED is an input argument if FACT = 'F'; otherwise, it is an */
	/*          output argument. */

	/*  R       (input or output) DOUBLE PRECISION array, dimension (N) */
	/*          The row scale factors for A.  If EQUED = 'R' or 'B', A is */
	/*          multiplied on the left by diag(R); if EQUED = 'N' or 'C', R */
	/*          is not accessed.  R is an input argument if FACT = 'F'; */
	/*          otherwise, R is an output argument.  If FACT = 'F' and */
	/*          EQUED = 'R' or 'B', each element of R must be positive. */

	/*  C       (input or output) DOUBLE PRECISION array, dimension (N) */
	/*          The column scale factors for A.  If EQUED = 'C' or 'B', A is */
	/*          multiplied on the right by diag(C); if EQUED = 'N' or 'R', C */
	/*          is not accessed.  C is an input argument if FACT = 'F'; */
	/*          otherwise, C is an output argument.  If FACT = 'F' and */
	/*          EQUED = 'C' or 'B', each element of C must be positive. */

	/*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS) */
	/*          On entry, the N-by-NRHS right hand side matrix B. */
	/*          On exit, */
	/*          if EQUED = 'N', B is not modified; */
	/*          if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by */
	/*          diag(R)*B; */
	/*          if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is */
	/*          overwritten by diag(C)*B. */

	/*  LDB     (input) INTEGER */
	/*          The leading dimension of the array B.  LDB >= max(1,N). */

	/*  X       (output) DOUBLE PRECISION array, dimension (LDX,NRHS) */
	/*          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X */
	/*          to the original system of equations.  Note that A and B are */
	/*          modified on exit if EQUED .ne. 'N', and the solution to the */
	/*          equilibrated system is inv(diag(C))*X if TRANS = 'N' and */
	/*          EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C' */
	/*          and EQUED = 'R' or 'B'. */

	/*  LDX     (input) INTEGER */
	/*          The leading dimension of the array X.  LDX >= max(1,N). */

	/*  RCOND   (output) DOUBLE PRECISION */
	/*          The estimate of the reciprocal condition number of the matrix */
	/*          A after equilibration (if done).  If RCOND is less than the */
	/*          machine precision (in particular, if RCOND = 0), the matrix */
	/*          is singular to working precision.  This condition is */
	/*          indicated by a return code of INFO > 0. */

	/*  FERR    (output) DOUBLE PRECISION array, dimension (NRHS) */
	/*          The estimated forward error bound for each solution vector */
	/*          X(j) (the j-th column of the solution matrix X). */
	/*          If XTRUE is the true solution corresponding to X(j), FERR(j) */
	/*          is an estimated upper bound for the magnitude of the largest */
	/*          element in (X(j) - XTRUE) divided by the magnitude of the */
	/*          largest element in X(j).  The estimate is as reliable as */
	/*          the estimate for RCOND, and is almost always a slight */
	/*          overestimate of the true error. */

	/*  BERR    (output) DOUBLE PRECISION array, dimension (NRHS) */
	/*          The componentwise relative backward error of each solution */
	/*          vector X(j) (i.e., the smallest relative change in */
	/*          any element of A or B that makes X(j) an exact solution). */

	/*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (4*N) */
	/*          On exit, WORK(1) contains the reciprocal pivot growth */
	/*          factor norm(A)/norm(U). The "max absolute element" norm is */
	/*          used. If WORK(1) is much less than 1, then the stability */
	/*          of the LU factorization of the (equilibrated) matrix A */
	/*          could be poor. This also means that the solution X, condition */
	/*          estimator RCOND, and forward error bound FERR could be */
	/*          unreliable. If factorization fails with 0<INFO<=N, then */
	/*          WORK(1) contains the reciprocal pivot growth factor for the */
	/*          leading INFO columns of A. */

	/*  IWORK   (workspace) INTEGER array, dimension (N) */

	/*  INFO    (output) INTEGER */
	/*          = 0:  successful exit */
	/*          < 0:  if INFO = -i, the i-th argument had an illegal value */
	/*          > 0:  if INFO = i, and i is */
	/*                <= N:  U(i,i) is exactly zero.  The factorization has */
	/*                       been completed, but the factor U is exactly */
	/*                       singular, so the solution and error bounds */
	/*                       could not be computed. RCOND = 0 is returned. */
	/*                = N+1: U is nonsingular, but RCOND is less than machine */
	/*                       precision, meaning that the matrix is singular */
	/*                       to working precision.  Nevertheless, the */
	/*                       solution and error bounds are computed because */
	/*                       there are a number of situations where the */
	/*                       computed solution can be more accurate than the */
	/*                       value of RCOND would suggest. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	af_dim1 = *ldaf;
	af_offset = 1 + af_dim1;
	af -= af_offset;
	--ipiv;
	--r__;
	--c__;
	b_dim1 = *ldb;
	b_offset = 1 + b_dim1;
	b -= b_offset;
	x_dim1 = *ldx;
	x_offset = 1 + x_dim1;
	x -= x_offset;
	--ferr;
	--berr;
	--work;
	--iwork;

	/* Function Body */
	*info = 0;
	nofact = lsame_ (fact, "N", (ftnlen) 1, (ftnlen) 1);
	equil = lsame_ (fact, "E", (ftnlen) 1, (ftnlen) 1);
	notran = lsame_ (trans, "N", (ftnlen) 1, (ftnlen) 1);
	if (nofact || equil)
	{
		*(unsigned char *) equed = 'N';
		rowequ = FALSE_;
		colequ = FALSE_;
	}
	else
	{
		rowequ = lsame_ (equed, "R", (ftnlen) 1, (ftnlen) 1) || lsame_ (equed, "B", (ftnlen) 1, (ftnlen) 1);
		colequ = lsame_ (equed, "C", (ftnlen) 1, (ftnlen) 1) || lsame_ (equed, "B", (ftnlen) 1, (ftnlen) 1);
		smlnum = dlamch_ ("Safe minimum", (ftnlen) 12);
		bignum = 1. / smlnum;
	}

	/*     Test the input parameters. */

	if (!nofact && !equil && !lsame_ (fact, "F", (ftnlen) 1, (ftnlen) 1))
	{
		*info = -1;
	}
	else if (!notran && !lsame_ (trans, "T", (ftnlen) 1, (ftnlen) 1) && !lsame_ (trans, "C", (ftnlen) 1, (ftnlen) 1))
	{
		*info = -2;
	}
	else if (*n < 0)
	{
		*info = -3;
	}
	else if (*nrhs < 0)
	{
		*info = -4;
	}
	else if (*lda < max (1, *n))
	{
		*info = -6;
	}
	else if (*ldaf < max (1, *n))
	{
		*info = -8;
	}
	else if (lsame_ (fact, "F", (ftnlen) 1, (ftnlen) 1) && !(rowequ || colequ || lsame_ (equed, "N", (ftnlen) 1, (ftnlen) 1)))
	{
		*info = -10;
	}
	else
	{
		if (rowequ)
		{
			rcmin = bignum;
			rcmax = 0.;
			i__1 = *n;
			for (j = 1; j <= i__1; ++j)
			{
				/* Computing MIN */
				d__1 = rcmin, d__2 = r__[j];
				rcmin = min (d__1, d__2);
				/* Computing MAX */
				d__1 = rcmax, d__2 = r__[j];
				rcmax = max (d__1, d__2);
				/* L10: */
			}
			if (rcmin <= 0.)
			{
				*info = -11;
			}
			else if (*n > 0)
			{
				rowcnd = max (rcmin, smlnum) / min (rcmax, bignum);
			}
			else
			{
				rowcnd = 1.;
			}
		}
		if (colequ && *info == 0)
		{
			rcmin = bignum;
			rcmax = 0.;
			i__1 = *n;
			for (j = 1; j <= i__1; ++j)
			{
				/* Computing MIN */
				d__1 = rcmin, d__2 = c__[j];
				rcmin = min (d__1, d__2);
				/* Computing MAX */
				d__1 = rcmax, d__2 = c__[j];
				rcmax = max (d__1, d__2);
				/* L20: */
			}
			if (rcmin <= 0.)
			{
				*info = -12;
			}
			else if (*n > 0)
			{
				colcnd = max (rcmin, smlnum) / min (rcmax, bignum);
			}
			else
			{
				colcnd = 1.;
			}
		}
		if (*info == 0)
		{
			if (*ldb < max (1, *n))
			{
				*info = -14;
			}
			else if (*ldx < max (1, *n))
			{
				*info = -16;
			}
		}
	}

	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DGESVX", &i__1, (ftnlen) 6);
		return 0;
	}

	if (equil)
	{

		/*        Compute row and column scalings to equilibrate the matrix A. */

		dgeequ_ (n, n, &a[a_offset], lda, &r__[1], &c__[1], &rowcnd, &colcnd, &amax, &infequ);
		if (infequ == 0)
		{

			/*           Equilibrate the matrix. */

			dlaqge_ (n, n, &a[a_offset], lda, &r__[1], &c__[1], &rowcnd, &colcnd, &amax, equed, (ftnlen) 1);
			rowequ = lsame_ (equed, "R", (ftnlen) 1, (ftnlen) 1) || lsame_ (equed, "B", (ftnlen) 1, (ftnlen) 1);
			colequ = lsame_ (equed, "C", (ftnlen) 1, (ftnlen) 1) || lsame_ (equed, "B", (ftnlen) 1, (ftnlen) 1);
		}
	}

	/*     Scale the right hand side. */

	if (notran)
	{
		if (rowequ)
		{
			i__1 = *nrhs;
			for (j = 1; j <= i__1; ++j)
			{
				i__2 = *n;
				for (i__ = 1; i__ <= i__2; ++i__)
				{
					b[i__ + j * b_dim1] = r__[i__] * b[i__ + j * b_dim1];
					/* L30: */
				}
				/* L40: */
			}
		}
	}
	else if (colequ)
	{
		i__1 = *nrhs;
		for (j = 1; j <= i__1; ++j)
		{
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				b[i__ + j * b_dim1] = c__[i__] * b[i__ + j * b_dim1];
				/* L50: */
			}
			/* L60: */
		}
	}

	if (nofact || equil)
	{

		/*        Compute the LU factorization of A. */

		dlacpy_ ("Full", n, n, &a[a_offset], lda, &af[af_offset], ldaf, (ftnlen) 4);
		dgetrf_ (n, n, &af[af_offset], ldaf, &ipiv[1], info);

		/*        Return if INFO is non-zero. */

		if (*info != 0)
		{
			if (*info > 0)
			{

				/*              Compute the reciprocal pivot growth factor of the */
				/*              leading rank-deficient INFO columns of A. */

				rpvgrw = dlantr_ ("M", "U", "N", info, info, &af[af_offset], ldaf, &work[1], (ftnlen) 1, (ftnlen) 1, (ftnlen) 1);
				if (rpvgrw == 0.)
				{
					rpvgrw = 1.;
				}
				else
				{
					rpvgrw = dlange_ ("M", n, info, &a[a_offset], lda, &work[1], (ftnlen) 1) / rpvgrw;
				}
				work[1] = rpvgrw;
				*rcond = 0.;
			}
			return 0;
		}
	}

	/*     Compute the norm of the matrix A and the */
	/*     reciprocal pivot growth factor RPVGRW. */

	if (notran)
	{
		*(unsigned char *) norm = '1';
	}
	else
	{
		*(unsigned char *) norm = 'I';
	}
	anorm = dlange_ (norm, n, n, &a[a_offset], lda, &work[1], (ftnlen) 1);
	rpvgrw = dlantr_ ("M", "U", "N", n, n, &af[af_offset], ldaf, &work[1], (ftnlen) 1, (ftnlen) 1, (ftnlen) 1);
	if (rpvgrw == 0.)
	{
		rpvgrw = 1.;
	}
	else
	{
		rpvgrw = dlange_ ("M", n, n, &a[a_offset], lda, &work[1], (ftnlen) 1) / rpvgrw;
	}

	/*     Compute the reciprocal of the condition number of A. */

	dgecon_ (norm, n, &af[af_offset], ldaf, &anorm, rcond, &work[1], &iwork[1], info, (ftnlen) 1);

	/*     Set INFO = N+1 if the matrix is singular to working precision. */

	if (*rcond < dlamch_ ("Epsilon", (ftnlen) 7))
	{
		*info = *n + 1;
	}

	/*     Compute the solution matrix X. */

	dlacpy_ ("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen) 4);
	dgetrs_ (trans, n, nrhs, &af[af_offset], ldaf, &ipiv[1], &x[x_offset], ldx, info, (ftnlen) 1);

	/*     Use iterative refinement to improve the computed solution and */
	/*     compute error bounds and backward error estimates for it. */

	dgerfs_ (trans, n, nrhs, &a[a_offset], lda, &af[af_offset], ldaf, &ipiv[1],
				&b[b_offset], ldb, &x[x_offset], ldx, &ferr[1], &berr[1], &work[1], &iwork[1], info, (ftnlen) 1);

	/*     Transform the solution matrix X to a solution of the original */
	/*     system. */

	if (notran)
	{
		if (colequ)
		{
			i__1 = *nrhs;
			for (j = 1; j <= i__1; ++j)
			{
				i__2 = *n;
				for (i__ = 1; i__ <= i__2; ++i__)
				{
					x[i__ + j * x_dim1] = c__[i__] * x[i__ + j * x_dim1];
					/* L70: */
				}
				/* L80: */
			}
			i__1 = *nrhs;
			for (j = 1; j <= i__1; ++j)
			{
				ferr[j] /= colcnd;
				/* L90: */
			}
		}
	}
	else if (rowequ)
	{
		i__1 = *nrhs;
		for (j = 1; j <= i__1; ++j)
		{
			i__2 = *n;
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				x[i__ + j * x_dim1] = r__[i__] * x[i__ + j * x_dim1];
				/* L100: */
			}
			/* L110: */
		}
		i__1 = *nrhs;
		for (j = 1; j <= i__1; ++j)
		{
			ferr[j] /= rowcnd;
			/* L120: */
		}
	}

	work[1] = rpvgrw;
	return 0;

	/*     End of DGESVX */

}	/* dgesvx_ */

static int
dgetf2_ (integer * m, integer * n, doublereal * a, integer * lda, integer * ipiv, integer * info)
{
	/* System generated locals */
	integer a_dim1, a_offset, i__1, i__2, i__3;
	doublereal d__1;

	/* Local variables */
	static integer j, jp;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     June 30, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DGETF2 computes an LU factorization of a general m-by-n matrix A */
	/*  using partial pivoting with row interchanges. */

	/*  The factorization has the form */
	/*     A = P * L * U */
	/*  where P is a permutation matrix, L is lower triangular with unit */
	/*  diagonal elements (lower trapezoidal if m > n), and U is upper */
	/*  triangular (upper trapezoidal if m < n). */

	/*  This is the right-looking Level 2 BLAS version of the algorithm. */

	/*  Arguments */
	/*  ========= */

	/*  M       (input) INTEGER */
	/*          The number of rows of the matrix A.  M >= 0. */

	/*  N       (input) INTEGER */
	/*          The number of columns of the matrix A.  N >= 0. */

	/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          On entry, the m by n matrix to be factored. */
	/*          On exit, the factors L and U from the factorization */
	/*          A = P*L*U; the unit diagonal elements of L are not stored. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(1,M). */

	/*  IPIV    (output) INTEGER array, dimension (min(M,N)) */
	/*          The pivot indices; for 1 <= i <= min(M,N), row i of the */
	/*          matrix was interchanged with row IPIV(i). */

	/*  INFO    (output) INTEGER */
	/*          = 0: successful exit */
	/*          < 0: if INFO = -k, the k-th argument had an illegal value */
	/*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization */
	/*               has been completed, but the factor U is exactly */
	/*               singular, and division by zero will occur if it is used */
	/*               to solve a system of equations. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Test the input parameters. */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--ipiv;

	/* Function Body */
	*info = 0;
	if (*m < 0)
	{
		*info = -1;
	}
	else if (*n < 0)
	{
		*info = -2;
	}
	else if (*lda < max (1, *m))
	{
		*info = -4;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DGETF2", &i__1, (ftnlen) 6);
		return 0;
	}

	/*     Quick return if possible */

	if (*m == 0 || *n == 0)
	{
		return 0;
	}

	i__1 = min (*m, *n);
	for (j = 1; j <= i__1; ++j)
	{

		/*        Find pivot and test for singularity. */

		i__2 = *m - j + 1;
		jp = j - 1 + idamax_ (&i__2, &a[j + j * a_dim1], &c__1);
		ipiv[j] = jp;
		if (a[jp + j * a_dim1] != 0.)
		{

			/*           Apply the interchange to columns 1:N. */

			if (jp != j)
			{
				dswap_ (n, &a[j + a_dim1], lda, &a[jp + a_dim1], lda);
			}

			/*           Compute elements J+1:M of J-th column. */

			if (j < *m)
			{
				i__2 = *m - j;
				d__1 = 1. / a[j + j * a_dim1];
				dscal_ (&i__2, &d__1, &a[j + 1 + j * a_dim1], &c__1);
			}

		}
		else if (*info == 0)
		{

			*info = j;
		}

		if (j < min (*m, *n))
		{

			/*           Update trailing submatrix. */

			i__2 = *m - j;
			i__3 = *n - j;
			dger_ (&i__2, &i__3, &c_b347, &a[j + 1 + j * a_dim1], &c__1, &a[j + (j + 1) * a_dim1], lda, &a[j + 1 + (j + 1) * a_dim1], lda);
		}
		/* L10: */
	}
	return 0;

	/*     End of DGETF2 */

}	/* dgetf2_ */

static int
dgetrf_ (integer * m, integer * n, doublereal * a, integer * lda, integer * ipiv, integer * info)
{
	/* System generated locals */
	integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;

	/* Local variables */
	static integer i__, j, jb, nb;
	static integer iinfo;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     March 31, 1993 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DGETRF computes an LU factorization of a general M-by-N matrix A */
	/*  using partial pivoting with row interchanges. */

	/*  The factorization has the form */
	/*     A = P * L * U */
	/*  where P is a permutation matrix, L is lower triangular with unit */
	/*  diagonal elements (lower trapezoidal if m > n), and U is upper */
	/*  triangular (upper trapezoidal if m < n). */

	/*  This is the right-looking Level 3 BLAS version of the algorithm. */

	/*  Arguments */
	/*  ========= */

	/*  M       (input) INTEGER */
	/*          The number of rows of the matrix A.  M >= 0. */

	/*  N       (input) INTEGER */
	/*          The number of columns of the matrix A.  N >= 0. */

	/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          On entry, the M-by-N matrix to be factored. */
	/*          On exit, the factors L and U from the factorization */
	/*          A = P*L*U; the unit diagonal elements of L are not stored. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(1,M). */

	/*  IPIV    (output) INTEGER array, dimension (min(M,N)) */
	/*          The pivot indices; for 1 <= i <= min(M,N), row i of the */
	/*          matrix was interchanged with row IPIV(i). */

	/*  INFO    (output) INTEGER */
	/*          = 0:  successful exit */
	/*          < 0:  if INFO = -i, the i-th argument had an illegal value */
	/*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization */
	/*                has been completed, but the factor U is exactly */
	/*                singular, and division by zero will occur if it is used */
	/*                to solve a system of equations. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Test the input parameters. */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--ipiv;

	/* Function Body */
	*info = 0;
	if (*m < 0)
	{
		*info = -1;
	}
	else if (*n < 0)
	{
		*info = -2;
	}
	else if (*lda < max (1, *m))
	{
		*info = -4;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DGETRF", &i__1, (ftnlen) 6);
		return 0;
	}

	/*     Quick return if possible */

	if (*m == 0 || *n == 0)
	{
		return 0;
	}

	/*     Determine the block size for this environment. */

	nb = ilaenv_ (&c__1, "DGETRF", " ", m, n, &c_n1, &c_n1, (ftnlen) 6, (ftnlen) 1);
	if (nb <= 1 || nb >= min (*m, *n))
	{

		/*        Use unblocked code. */

		dgetf2_ (m, n, &a[a_offset], lda, &ipiv[1], info);
	}
	else
	{

		/*        Use blocked code. */

		i__1 = min (*m, *n);
		i__2 = nb;
		for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2)
		{
			/* Computing MIN */
			i__3 = min (*m, *n) - j + 1;
			jb = min (i__3, nb);

			/*           Factor diagonal and subdiagonal blocks and test for exact */
			/*           singularity. */

			i__3 = *m - j + 1;
			dgetf2_ (&i__3, &jb, &a[j + j * a_dim1], lda, &ipiv[j], &iinfo);

			/*           Adjust INFO and the pivot indices. */

			if (*info == 0 && iinfo > 0)
			{
				*info = iinfo + j - 1;
			}
			/* Computing MIN */
			i__4 = *m, i__5 = j + jb - 1;
			i__3 = min (i__4, i__5);
			for (i__ = j; i__ <= i__3; ++i__)
			{
				ipiv[i__] = j - 1 + ipiv[i__];
				/* L10: */
			}

			/*           Apply interchanges to columns 1:J-1. */

			i__3 = j - 1;
			i__4 = j + jb - 1;
			dlaswp_ (&i__3, &a[a_offset], lda, &j, &i__4, &ipiv[1], &c__1);

			if (j + jb <= *n)
			{

				/*              Apply interchanges to columns J+JB:N. */

				i__3 = *n - j - jb + 1;
				i__4 = j + jb - 1;
				dlaswp_ (&i__3, &a[(j + jb) * a_dim1 + 1], lda, &j, &i__4, &ipiv[1], &c__1);

				/*              Compute block row of U. */

				i__3 = *n - j - jb + 1;
				dtrsm_ ("Left", "Lower", "No transpose", "Unit", &jb, &i__3, &c_b348, &a[j + j * a_dim1], lda, &a[j + (j + jb) *
																																				  a_dim1], lda, (ftnlen) 4, (ftnlen) 5, (ftnlen) 12,
						  (ftnlen) 4);
				if (j + jb <= *m)
				{

					/*                 Update trailing submatrix. */

					i__3 = *m - j - jb + 1;
					i__4 = *n - j - jb + 1;
					dgemm_ ("No transpose", "No transpose", &i__3, &i__4, &jb,
							  &c_b347, &a[j + jb + j * a_dim1], lda, &a[j + (j
																							 + jb) * a_dim1], lda, &c_b348, &a[j + jb + (j +
																																						jb) * a_dim1], lda, (ftnlen) 12, (ftnlen) 12);
				}
			}
			/* L20: */
		}
	}
	return 0;

	/*     End of DGETRF */

}	/* dgetrf_ */

static int
dgetrs_ (char *trans, integer * n, integer * nrhs,
			doublereal * a, integer * lda, integer * ipiv, doublereal * b, integer * ldb, integer * info, ftnlen trans_len)
{
	/* System generated locals */
	integer a_dim1, a_offset, b_dim1, b_offset, i__1;

	/* Local variables */
	static logical notran;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     March 31, 1993 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DGETRS solves a system of linear equations */
	/*     A * X = B  or  A' * X = B */
	/*  with a general N-by-N matrix A using the LU factorization computed */
	/*  by DGETRF. */

	/*  Arguments */
	/*  ========= */

	/*  TRANS   (input) CHARACTER*1 */
	/*          Specifies the form of the system of equations: */
	/*          = 'N':  A * X = B  (No transpose) */
	/*          = 'T':  A'* X = B  (Transpose) */
	/*          = 'C':  A'* X = B  (Conjugate transpose = Transpose) */

	/*  N       (input) INTEGER */
	/*          The order of the matrix A.  N >= 0. */

	/*  NRHS    (input) INTEGER */
	/*          The number of right hand sides, i.e., the number of columns */
	/*          of the matrix B.  NRHS >= 0. */

	/*  A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          The factors L and U from the factorization A = P*L*U */
	/*          as computed by DGETRF. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(1,N). */

	/*  IPIV    (input) INTEGER array, dimension (N) */
	/*          The pivot indices from DGETRF; for 1<=i<=N, row i of the */
	/*          matrix was interchanged with row IPIV(i). */

	/*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS) */
	/*          On entry, the right hand side matrix B. */
	/*          On exit, the solution matrix X. */

	/*  LDB     (input) INTEGER */
	/*          The leading dimension of the array B.  LDB >= max(1,N). */

	/*  INFO    (output) INTEGER */
	/*          = 0:  successful exit */
	/*          < 0:  if INFO = -i, the i-th argument had an illegal value */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Test the input parameters. */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--ipiv;
	b_dim1 = *ldb;
	b_offset = 1 + b_dim1;
	b -= b_offset;

	/* Function Body */
	*info = 0;
	notran = lsame_ (trans, "N", (ftnlen) 1, (ftnlen) 1);
	if (!notran && !lsame_ (trans, "T", (ftnlen) 1, (ftnlen) 1) && !lsame_ (trans, "C", (ftnlen) 1, (ftnlen) 1))
	{
		*info = -1;
	}
	else if (*n < 0)
	{
		*info = -2;
	}
	else if (*nrhs < 0)
	{
		*info = -3;
	}
	else if (*lda < max (1, *n))
	{
		*info = -5;
	}
	else if (*ldb < max (1, *n))
	{
		*info = -8;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DGETRS", &i__1, (ftnlen) 6);
		return 0;
	}

	/*     Quick return if possible */

	if (*n == 0 || *nrhs == 0)
	{
		return 0;
	}

	if (notran)
	{

		/*        Solve A * X = B. */

		/*        Apply row interchanges to the right hand sides. */

		dlaswp_ (nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c__1);

		/*        Solve L*X = B, overwriting B with X. */

		dtrsm_ ("Left", "Lower", "No transpose", "Unit", n, nrhs, &c_b348, &a[a_offset], lda, &b[b_offset], ldb, (ftnlen) 4, (ftnlen) 5, (ftnlen) 12, (ftnlen) 4);

		/*        Solve U*X = B, overwriting B with X. */

		dtrsm_ ("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &c_b348,
				  &a[a_offset], lda, &b[b_offset], ldb, (ftnlen) 4, (ftnlen) 5, (ftnlen) 12, (ftnlen) 8);
	}
	else
	{

		/*        Solve A' * X = B. */

		/*        Solve U'*X = B, overwriting B with X. */

		dtrsm_ ("Left", "Upper", "Transpose", "Non-unit", n, nrhs, &c_b348, &a[a_offset], lda, &b[b_offset], ldb, (ftnlen) 4, (ftnlen) 5, (ftnlen) 9, (ftnlen) 8);

		/*        Solve L'*X = B, overwriting B with X. */

		dtrsm_ ("Left", "Lower", "Transpose", "Unit", n, nrhs, &c_b348, &a[a_offset], lda, &b[b_offset], ldb, (ftnlen) 4, (ftnlen) 5, (ftnlen) 9, (ftnlen) 4);

		/*        Apply row interchanges to the solution vectors. */

		dlaswp_ (nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c_n1);
	}

	return 0;

	/*     End of DGETRS */

}	/* dgetrs_ */

static int
dhseqr_ (char *job, char *compz, integer * n, integer * ilo,
			integer * ihi, doublereal * h__, integer * ldh, doublereal * wr,
			doublereal * wi, doublereal * z__, integer * ldz, doublereal * work, integer * lwork, integer * info, ftnlen job_len, ftnlen compz_len)
{
	/* System generated locals */
	address a__1[2];
	integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3[2], i__4, i__5;
	doublereal d__1, d__2;
	char ch__1[2];

	/* Builtin functions */
	/* Subroutine */ int s_cat (char *, char **, integer *, integer *, ftnlen);

	/* Local variables */
	static integer i__, j, k, l;
	static doublereal s[225] /* was [15][15] */ , v[16];
	static integer i1, i2, ii, nh, nr, ns, nv;
	static doublereal vv[16];
	static integer itn;
	static doublereal tau;
	static integer its;
	static doublereal ulp, tst1;
	static integer maxb;
	static doublereal absw;
	static integer ierr;
	static doublereal unfl, temp, ovfl;
	static integer itemp;
	static logical initz, wantt, wantz;
	static doublereal smlnum;
	static logical lquery;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     June 30, 1999 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DHSEQR computes the eigenvalues of a real upper Hessenberg matrix H */
	/*  and, optionally, the matrices T and Z from the Schur decomposition */
	/*  H = Z T Z**T, where T is an upper quasi-triangular matrix (the Schur */
	/*  form), and Z is the orthogonal matrix of Schur vectors. */

	/*  Optionally Z may be postmultiplied into an input orthogonal matrix Q, */
	/*  so that this routine can give the Schur factorization of a matrix A */
	/*  which has been reduced to the Hessenberg form H by the orthogonal */
	/*  matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T. */

	/*  Arguments */
	/*  ========= */

	/*  JOB     (input) CHARACTER*1 */
	/*          = 'E':  compute eigenvalues only; */
	/*          = 'S':  compute eigenvalues and the Schur form T. */

	/*  COMPZ   (input) CHARACTER*1 */
	/*          = 'N':  no Schur vectors are computed; */
	/*          = 'I':  Z is initialized to the unit matrix and the matrix Z */
	/*                  of Schur vectors of H is returned; */
	/*          = 'V':  Z must contain an orthogonal matrix Q on entry, and */
	/*                  the product Q*Z is returned. */

	/*  N       (input) INTEGER */
	/*          The order of the matrix H.  N >= 0. */

	/*  ILO     (input) INTEGER */
	/*  IHI     (input) INTEGER */
	/*          It is assumed that H is already upper triangular in rows */
	/*          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally */
	/*          set by a previous call to DGEBAL, and then passed to SGEHRD */
	/*          when the matrix output by DGEBAL is reduced to Hessenberg */
	/*          form. Otherwise ILO and IHI should be set to 1 and N */
	/*          respectively. */
	/*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */

	/*  H       (input/output) DOUBLE PRECISION array, dimension (LDH,N) */
	/*          On entry, the upper Hessenberg matrix H. */
	/*          On exit, if JOB = 'S', H contains the upper quasi-triangular */
	/*          matrix T from the Schur decomposition (the Schur form); */
	/*          2-by-2 diagonal blocks (corresponding to complex conjugate */
	/*          pairs of eigenvalues) are returned in standard form, with */
	/*          H(i,i) = H(i+1,i+1) and H(i+1,i)*H(i,i+1) < 0. If JOB = 'E', */
	/*          the contents of H are unspecified on exit. */

	/*  LDH     (input) INTEGER */
	/*          The leading dimension of the array H. LDH >= max(1,N). */

	/*  WR      (output) DOUBLE PRECISION array, dimension (N) */
	/*  WI      (output) DOUBLE PRECISION array, dimension (N) */
	/*          The real and imaginary parts, respectively, of the computed */
	/*          eigenvalues. If two eigenvalues are computed as a complex */
	/*          conjugate pair, they are stored in consecutive elements of */
	/*          WR and WI, say the i-th and (i+1)th, with WI(i) > 0 and */
	/*          WI(i+1) < 0. If JOB = 'S', the eigenvalues are stored in the */
	/*          same order as on the diagonal of the Schur form returned in */
	/*          H, with WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2 */
	/*          diagonal block, WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and */
	/*          WI(i+1) = -WI(i). */

	/*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N) */
	/*          If COMPZ = 'N': Z is not referenced. */
	/*          If COMPZ = 'I': on entry, Z need not be set, and on exit, Z */
	/*          contains the orthogonal matrix Z of the Schur vectors of H. */
	/*          If COMPZ = 'V': on entry Z must contain an N-by-N matrix Q, */
	/*          which is assumed to be equal to the unit matrix except for */
	/*          the submatrix Z(ILO:IHI,ILO:IHI); on exit Z contains Q*Z. */
	/*          Normally Q is the orthogonal matrix generated by DORGHR after */
	/*          the call to DGEHRD which formed the Hessenberg matrix H. */

	/*  LDZ     (input) INTEGER */
	/*          The leading dimension of the array Z. */
	/*          LDZ >= max(1,N) if COMPZ = 'I' or 'V'; LDZ >= 1 otherwise. */

	/*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) */
	/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

	/*  LWORK   (input) INTEGER */
	/*          The dimension of the array WORK.  LWORK >= max(1,N). */

	/*          If LWORK = -1, then a workspace query is assumed; the routine */
	/*          only calculates the optimal size of the WORK array, returns */
	/*          this value as the first entry of the WORK array, and no error */
	/*          message related to LWORK is issued by XERBLA. */

	/*  INFO    (output) INTEGER */
	/*          = 0:  successful exit */
	/*          < 0:  if INFO = -i, the i-th argument had an illegal value */
	/*          > 0:  if INFO = i, DHSEQR failed to compute all of the */
	/*                eigenvalues in a total of 30*(IHI-ILO+1) iterations; */
	/*                elements 1:ilo-1 and i+1:n of WR and WI contain those */
	/*                eigenvalues which have been successfully computed. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Local Arrays .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Decode and test the input parameters */

	/* Parameter adjustments */
	h_dim1 = *ldh;
	h_offset = 1 + h_dim1;
	h__ -= h_offset;
	--wr;
	--wi;
	z_dim1 = *ldz;
	z_offset = 1 + z_dim1;
	z__ -= z_offset;
	--work;

	/* Function Body */
	wantt = lsame_ (job, "S", (ftnlen) 1, (ftnlen) 1);
	initz = lsame_ (compz, "I", (ftnlen) 1, (ftnlen) 1);
	wantz = initz || lsame_ (compz, "V", (ftnlen) 1, (ftnlen) 1);

	*info = 0;
	work[1] = (doublereal) max (1, *n);
	lquery = *lwork == -1;
	if (!lsame_ (job, "E", (ftnlen) 1, (ftnlen) 1) && !wantt)
	{
		*info = -1;
	}
	else if (!lsame_ (compz, "N", (ftnlen) 1, (ftnlen) 1) && !wantz)
	{
		*info = -2;
	}
	else if (*n < 0)
	{
		*info = -3;
	}
	else if (*ilo < 1 || *ilo > max (1, *n))
	{
		*info = -4;
	}
	else if (*ihi < min (*ilo, *n) || *ihi > *n)
	{
		*info = -5;
	}
	else if (*ldh < max (1, *n))
	{
		*info = -7;
	}
	else if (*ldz < 1 || wantz && *ldz < max (1, *n))
	{
		*info = -11;
	}
	else if (*lwork < max (1, *n) && !lquery)
	{
		*info = -13;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DHSEQR", &i__1, (ftnlen) 6);
		return 0;
	}
	else if (lquery)
	{
		return 0;
	}

	/*     Initialize Z, if necessary */

	if (initz)
	{
		dlaset_ ("Full", n, n, &c_b507, &c_b348, &z__[z_offset], ldz, (ftnlen) 4);
	}

	/*     Store the eigenvalues isolated by DGEBAL. */

	i__1 = *ilo - 1;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		wr[i__] = h__[i__ + i__ * h_dim1];
		wi[i__] = 0.;
		/* L10: */
	}
	i__1 = *n;
	for (i__ = *ihi + 1; i__ <= i__1; ++i__)
	{
		wr[i__] = h__[i__ + i__ * h_dim1];
		wi[i__] = 0.;
		/* L20: */
	}

	/*     Quick return if possible. */

	if (*n == 0)
	{
		return 0;
	}
	if (*ilo == *ihi)
	{
		wr[*ilo] = h__[*ilo + *ilo * h_dim1];
		wi[*ilo] = 0.;
		return 0;
	}

	/*     Set rows and columns ILO to IHI to zero below the first */
	/*     subdiagonal. */

	i__1 = *ihi - 2;
	for (j = *ilo; j <= i__1; ++j)
	{
		i__2 = *n;
		for (i__ = j + 2; i__ <= i__2; ++i__)
		{
			h__[i__ + j * h_dim1] = 0.;
			/* L30: */
		}
		/* L40: */
	}
	nh = *ihi - *ilo + 1;

	/*     Determine the order of the multi-shift QR algorithm to be used. */

	/* Writing concatenation */
	i__3[0] = 1, a__1[0] = job;
	i__3[1] = 1, a__1[1] = compz;
	s_cat (ch__1, a__1, i__3, &c__2, (ftnlen) 2);
	ns = ilaenv_ (&c__4, "DHSEQR", ch__1, n, ilo, ihi, &c_n1, (ftnlen) 6, (ftnlen) 2);
	/* Writing concatenation */
	i__3[0] = 1, a__1[0] = job;
	i__3[1] = 1, a__1[1] = compz;
	s_cat (ch__1, a__1, i__3, &c__2, (ftnlen) 2);
	maxb = ilaenv_ (&c__8, "DHSEQR", ch__1, n, ilo, ihi, &c_n1, (ftnlen) 6, (ftnlen) 2);
	if (ns <= 2 || ns > nh || maxb >= nh)
	{

		/*        Use the standard double-shift algorithm */

		dlahqr_ (&wantt, &wantz, n, ilo, ihi, &h__[h_offset], ldh, &wr[1], &wi[1], ilo, ihi, &z__[z_offset], ldz, info);
		return 0;
	}
	maxb = max (3, maxb);
	/* Computing MIN */
	i__1 = min (ns, maxb);
	ns = min (i__1, 15);

	/*     Now 2 < NS <= MAXB < NH. */

	/*     Set machine-dependent constants for the stopping criterion. */
	/*     If norm(H) <= sqrt(OVFL), overflow should not occur. */

	unfl = dlamch_ ("Safe minimum", (ftnlen) 12);
	ovfl = 1. / unfl;
	dlabad_ (&unfl, &ovfl);
	ulp = dlamch_ ("Precision", (ftnlen) 9);
	smlnum = unfl * (nh / ulp);

	/*     I1 and I2 are the indices of the first row and last column of H */
	/*     to which transformations must be applied. If eigenvalues only are */
	/*     being computed, I1 and I2 are set inside the main loop. */

	if (wantt)
	{
		i1 = 1;
		i2 = *n;
	}

	/*     ITN is the total number of multiple-shift QR iterations allowed. */

	itn = nh * 30;

	/*     The main loop begins here. I is the loop index and decreases from */
	/*     IHI to ILO in steps of at most MAXB. Each iteration of the loop */
	/*     works with the active submatrix in rows and columns L to I. */
	/*     Eigenvalues I+1 to IHI have already converged. Either L = ILO or */
	/*     H(L,L-1) is negligible so that the matrix splits. */

	i__ = *ihi;
 L50:
	l = *ilo;
	if (i__ < *ilo)
	{
		goto L170;
	}

	/*     Perform multiple-shift QR iterations on rows and columns ILO to I */
	/*     until a submatrix of order at most MAXB splits off at the bottom */
	/*     because a subdiagonal element has become negligible. */

	i__1 = itn;
	for (its = 0; its <= i__1; ++its)
	{

		/*        Look for a single small subdiagonal element. */

		i__2 = l + 1;
		for (k = i__; k >= i__2; --k)
		{
			tst1 = (d__1 = h__[k - 1 + (k - 1) * h_dim1], abs (d__1)) + (d__2 = h__[k + k * h_dim1], abs (d__2));
			if (tst1 == 0.)
			{
				i__4 = i__ - l + 1;
				tst1 = dlanhs_ ("1", &i__4, &h__[l + l * h_dim1], ldh, &work[1], (ftnlen) 1);
			}
			/* Computing MAX */
			d__2 = ulp * tst1;
			if ((d__1 = h__[k + (k - 1) * h_dim1], abs (d__1)) <= max (d__2, smlnum))
			{
				goto L70;
			}
			/* L60: */
		}
	 L70:
		l = k;
		if (l > *ilo)
		{

			/*           H(L,L-1) is negligible. */

			h__[l + (l - 1) * h_dim1] = 0.;
		}

		/*        Exit from loop if a submatrix of order <= MAXB has split off. */

		if (l >= i__ - maxb + 1)
		{
			goto L160;
		}

		/*        Now the active submatrix is in rows and columns L to I. If */
		/*        eigenvalues only are being computed, only the active submatrix */
		/*        need be transformed. */

		if (!wantt)
		{
			i1 = l;
			i2 = i__;
		}

		if (its == 20 || its == 30)
		{

			/*           Exceptional shifts. */

			i__2 = i__;
			for (ii = i__ - ns + 1; ii <= i__2; ++ii)
			{
				wr[ii] = ((d__1 = h__[ii + (ii - 1) * h_dim1], abs (d__1)) + (d__2 = h__[ii + ii * h_dim1], abs (d__2))) * 1.5;
				wi[ii] = 0.;
				/* L80: */
			}
		}
		else
		{

			/*           Use eigenvalues of trailing submatrix of order NS as shifts. */

			dlacpy_ ("Full", &ns, &ns, &h__[i__ - ns + 1 + (i__ - ns + 1) * h_dim1], ldh, s, &c__15, (ftnlen) 4);
			dlahqr_ (&c_false, &c_false, &ns, &c__1, &ns, s, &c__15, &wr[i__ - ns + 1], &wi[i__ - ns + 1], &c__1, &ns, &z__[z_offset], ldz, &ierr);
			if (ierr > 0)
			{

				/*              If DLAHQR failed to compute all NS eigenvalues, use the */
				/*              unconverged diagonal elements as the remaining shifts. */

				i__2 = ierr;
				for (ii = 1; ii <= i__2; ++ii)
				{
					wr[i__ - ns + ii] = s[ii + ii * 15 - 16];
					wi[i__ - ns + ii] = 0.;
					/* L90: */
				}
			}
		}

		/*        Form the first column of (G-w(1)) (G-w(2)) . . . (G-w(ns)) */
		/*        where G is the Hessenberg submatrix H(L:I,L:I) and w is */
		/*        the vector of shifts (stored in WR and WI). The result is */
		/*        stored in the local array V. */

		v[0] = 1.;
		i__2 = ns + 1;
		for (ii = 2; ii <= i__2; ++ii)
		{
			v[ii - 1] = 0.;
			/* L100: */
		}
		nv = 1;
		i__2 = i__;
		for (j = i__ - ns + 1; j <= i__2; ++j)
		{
			if (wi[j] >= 0.)
			{
				if (wi[j] == 0.)
				{

					/*                 real shift */

					i__4 = nv + 1;
					dcopy_ (&i__4, v, &c__1, vv, &c__1);
					i__4 = nv + 1;
					d__1 = -wr[j];
					dgemv_ ("No transpose", &i__4, &nv, &c_b348, &h__[l + l * h_dim1], ldh, vv, &c__1, &d__1, v, &c__1, (ftnlen) 12);
					++nv;
				}
				else if (wi[j] > 0.)
				{

					/*                 complex conjugate pair of shifts */

					i__4 = nv + 1;
					dcopy_ (&i__4, v, &c__1, vv, &c__1);
					i__4 = nv + 1;
					d__1 = wr[j] * -2.;
					dgemv_ ("No transpose", &i__4, &nv, &c_b348, &h__[l + l * h_dim1], ldh, v, &c__1, &d__1, vv, &c__1, (ftnlen) 12);
					i__4 = nv + 1;
					itemp = idamax_ (&i__4, vv, &c__1);
					/* Computing MAX */
					d__2 = (d__1 = vv[itemp - 1], abs (d__1));
					temp = 1. / max (d__2, smlnum);
					i__4 = nv + 1;
					dscal_ (&i__4, &temp, vv, &c__1);
					absw = dlapy2_ (&wr[j], &wi[j]);
					temp = temp * absw * absw;
					i__4 = nv + 2;
					i__5 = nv + 1;
					dgemv_ ("No transpose", &i__4, &i__5, &c_b348, &h__[l + l * h_dim1], ldh, vv, &c__1, &temp, v, &c__1, (ftnlen) 12);
					nv += 2;
				}

				/*              Scale V(1:NV) so that max(abs(V(i))) = 1. If V is zero, */
				/*              reset it to the unit vector. */

				itemp = idamax_ (&nv, v, &c__1);
				temp = (d__1 = v[itemp - 1], abs (d__1));
				if (temp == 0.)
				{
					v[0] = 1.;
					i__4 = nv;
					for (ii = 2; ii <= i__4; ++ii)
					{
						v[ii - 1] = 0.;
						/* L110: */
					}
				}
				else
				{
					temp = max (temp, smlnum);
					d__1 = 1. / temp;
					dscal_ (&nv, &d__1, v, &c__1);
				}
			}
			/* L120: */
		}

		/*        Multiple-shift QR step */

		i__2 = i__ - 1;
		for (k = l; k <= i__2; ++k)
		{

			/*           The first iteration of this loop determines a reflection G */
			/*           from the vector V and applies it from left and right to H, */
			/*           thus creating a nonzero bulge below the subdiagonal. */

			/*           Each subsequent iteration determines a reflection G to */
			/*           restore the Hessenberg form in the (K-1)th column, and thus */
			/*           chases the bulge one step toward the bottom of the active */
			/*           submatrix. NR is the order of G. */

			/* Computing MIN */
			i__4 = ns + 1, i__5 = i__ - k + 1;
			nr = min (i__4, i__5);
			if (k > l)
			{
				dcopy_ (&nr, &h__[k + (k - 1) * h_dim1], &c__1, v, &c__1);
			}
			dlarfg_ (&nr, v, &v[1], &c__1, &tau);
			if (k > l)
			{
				h__[k + (k - 1) * h_dim1] = v[0];
				i__4 = i__;
				for (ii = k + 1; ii <= i__4; ++ii)
				{
					h__[ii + (k - 1) * h_dim1] = 0.;
					/* L130: */
				}
			}
			v[0] = 1.;

			/*           Apply G from the left to transform the rows of the matrix in */
			/*           columns K to I2. */

			i__4 = i2 - k + 1;
			dlarfx_ ("Left", &nr, &i__4, v, &tau, &h__[k + k * h_dim1], ldh, &work[1], (ftnlen) 4);

			/*           Apply G from the right to transform the columns of the */
			/*           matrix in rows I1 to min(K+NR,I). */

			/* Computing MIN */
			i__5 = k + nr;
			i__4 = min (i__5, i__) - i1 + 1;
			dlarfx_ ("Right", &i__4, &nr, v, &tau, &h__[i1 + k * h_dim1], ldh, &work[1], (ftnlen) 5);

			if (wantz)
			{

				/*              Accumulate transformations in the matrix Z */

				dlarfx_ ("Right", &nh, &nr, v, &tau, &z__[*ilo + k * z_dim1], ldz, &work[1], (ftnlen) 5);
			}
			/* L140: */
		}

		/* L150: */
	}

	/*     Failure to converge in remaining number of iterations */

	*info = i__;
	return 0;

 L160:

	/*     A submatrix of order <= MAXB in rows and columns L to I has split */
	/*     off. Use the double-shift QR algorithm to handle it. */

	dlahqr_ (&wantt, &wantz, n, &l, &i__, &h__[h_offset], ldh, &wr[1], &wi[1], ilo, ihi, &z__[z_offset], ldz, info);
	if (*info > 0)
	{
		return 0;
	}

	/*     Decrement number of remaining iterations, and return to start of */
	/*     the main loop with a new value of I. */

	itn -= its;
	i__ = l - 1;
	goto L50;

 L170:
	work[1] = (doublereal) max (1, *n);
	return 0;

	/*     End of DHSEQR */

}	/* dhseqr_ */

static int
dlabad_ (doublereal * small, doublereal * large)
{
	/* Builtin functions */
	double d_lg10 (doublereal *), sqrt (doublereal);


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     October 31, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLABAD takes as input the values computed by DLAMCH for underflow and */
	/*  overflow, and returns the square root of each of these values if the */
	/*  log of LARGE is sufficiently large.  This subroutine is intended to */
	/*  identify machines with a large exponent range, such as the Crays, and */
	/*  redefine the underflow and overflow limits to be the square roots of */
	/*  the values computed by DLAMCH.  This subroutine is needed because */
	/*  DLAMCH does not compensate for poor arithmetic in the upper half of */
	/*  the exponent range, as is found on a Cray. */

	/*  Arguments */
	/*  ========= */

	/*  SMALL   (input/output) DOUBLE PRECISION */
	/*          On entry, the underflow threshold as computed by DLAMCH. */
	/*          On exit, if LOG10(LARGE) is sufficiently large, the square */
	/*          root of SMALL, otherwise unchanged. */

	/*  LARGE   (input/output) DOUBLE PRECISION */
	/*          On entry, the overflow threshold as computed by DLAMCH. */
	/*          On exit, if LOG10(LARGE) is sufficiently large, the square */
	/*          root of LARGE, otherwise unchanged. */

	/*  ===================================================================== */

	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     If it looks like we're on a Cray, take the square root of */
	/*     SMALL and LARGE to avoid overflow and underflow problems. */

	if (d_lg10 (large) > 2e3)
	{
		*small = sqrt (*small);
		*large = sqrt (*large);
	}

	return 0;

	/*     End of DLABAD */

}	/* dlabad_ */

static int
dlacon_ (integer * n, doublereal * v, doublereal * x, integer * isgn, doublereal * est, integer * kase)
{
	/* System generated locals */
	integer i__1;
	doublereal d__1;

	/* Builtin functions */
	double d_sign (doublereal *, doublereal *);
	integer i_dnnt (doublereal *);

	/* Local variables */
	static integer i__, j, iter;
	static doublereal temp;
	static integer jump;
	static integer jlast;
	static doublereal altsgn, estold;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     February 29, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLACON estimates the 1-norm of a square, real matrix A. */
	/*  Reverse communication is used for evaluating matrix-vector products. */

	/*  Arguments */
	/*  ========= */

	/*  N      (input) INTEGER */
	/*         The order of the matrix.  N >= 1. */

	/*  V      (workspace) DOUBLE PRECISION array, dimension (N) */
	/*         On the final return, V = A*W,  where  EST = norm(V)/norm(W) */
	/*         (W is not returned). */

	/*  X      (input/output) DOUBLE PRECISION array, dimension (N) */
	/*         On an intermediate return, X should be overwritten by */
	/*               A * X,   if KASE=1, */
	/*               A' * X,  if KASE=2, */
	/*         and DLACON must be re-called with all the other parameters */
	/*         unchanged. */

	/*  ISGN   (workspace) INTEGER array, dimension (N) */

	/*  EST    (output) DOUBLE PRECISION */
	/*         An estimate (a lower bound) for norm(A). */

	/*  KASE   (input/output) INTEGER */
	/*         On the initial call to DLACON, KASE should be 0. */
	/*         On an intermediate return, KASE will be 1 or 2, indicating */
	/*         whether X should be overwritten by A * X  or A' * X. */
	/*         On the final return from DLACON, KASE will again be 0. */

	/*  Further Details */
	/*  ======= ======= */

	/*  Contributed by Nick Higham, University of Manchester. */
	/*  Originally named SONEST, dated March 16, 1988. */

	/*  Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of */
	/*  a real or complex matrix, with applications to condition estimation", */
	/*  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Save statement .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/* Parameter adjustments */
	--isgn;
	--x;
	--v;

	/* Function Body */
	if (*kase == 0)
	{
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			x[i__] = 1. / (doublereal) (*n);
			/* L10: */
		}
		*kase = 1;
		jump = 1;
		return 0;
	}

	switch (jump)
	{
	case 1:
		goto L20;
	case 2:
		goto L40;
	case 3:
		goto L70;
	case 4:
		goto L110;
	case 5:
		goto L140;
	}

	/*     ................ ENTRY   (JUMP = 1) */
	/*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X. */

 L20:
	if (*n == 1)
	{
		v[1] = x[1];
		*est = abs (v[1]);
		/*        ... QUIT */
		goto L150;
	}
	*est = dasum_ (n, &x[1], &c__1);

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		x[i__] = d_sign (&c_b348, &x[i__]);
		isgn[i__] = i_dnnt (&x[i__]);
		/* L30: */
	}
	*kase = 2;
	jump = 2;
	return 0;

	/*     ................ ENTRY   (JUMP = 2) */
	/*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X. */

 L40:
	j = idamax_ (n, &x[1], &c__1);
	iter = 2;

	/*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */

 L50:
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		x[i__] = 0.;
		/* L60: */
	}
	x[j] = 1.;
	*kase = 1;
	jump = 3;
	return 0;

	/*     ................ ENTRY   (JUMP = 3) */
	/*     X HAS BEEN OVERWRITTEN BY A*X. */

 L70:
	dcopy_ (n, &x[1], &c__1, &v[1], &c__1);
	estold = *est;
	*est = dasum_ (n, &v[1], &c__1);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		d__1 = d_sign (&c_b348, &x[i__]);
		if (i_dnnt (&d__1) != isgn[i__])
		{
			goto L90;
		}
		/* L80: */
	}
	/*     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED. */
	goto L120;

 L90:
	/*     TEST FOR CYCLING. */
	if (*est <= estold)
	{
		goto L120;
	}

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		x[i__] = d_sign (&c_b348, &x[i__]);
		isgn[i__] = i_dnnt (&x[i__]);
		/* L100: */
	}
	*kase = 2;
	jump = 4;
	return 0;

	/*     ................ ENTRY   (JUMP = 4) */
	/*     X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X. */

 L110:
	jlast = j;
	j = idamax_ (n, &x[1], &c__1);
	if (x[jlast] != (d__1 = x[j], abs (d__1)) && iter < 5)
	{
		++iter;
		goto L50;
	}

	/*     ITERATION COMPLETE.  FINAL STAGE. */

 L120:
	altsgn = 1.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		x[i__] = altsgn * ((doublereal) (i__ - 1) / (doublereal) (*n - 1) + 1.);
		altsgn = -altsgn;
		/* L130: */
	}
	*kase = 1;
	jump = 5;
	return 0;

	/*     ................ ENTRY   (JUMP = 5) */
	/*     X HAS BEEN OVERWRITTEN BY A*X. */

 L140:
	temp = dasum_ (n, &x[1], &c__1) / (doublereal) (*n * 3) * 2.;
	if (temp > *est)
	{
		dcopy_ (n, &x[1], &c__1, &v[1], &c__1);
		*est = temp;
	}

 L150:
	*kase = 0;
	return 0;

	/*     End of DLACON */

}	/* dlacon_ */

static int
dlacpy_ (char *uplo, integer * m, integer * n, doublereal * a, integer * lda, doublereal * b, integer * ldb, ftnlen uplo_len)
{
	/* System generated locals */
	integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

	/* Local variables */
	static integer i__, j;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     February 29, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLACPY copies all or part of a two-dimensional matrix A to another */
	/*  matrix B. */

	/*  Arguments */
	/*  ========= */

	/*  UPLO    (input) CHARACTER*1 */
	/*          Specifies the part of the matrix A to be copied to B. */
	/*          = 'U':      Upper triangular part */
	/*          = 'L':      Lower triangular part */
	/*          Otherwise:  All of the matrix A */

	/*  M       (input) INTEGER */
	/*          The number of rows of the matrix A.  M >= 0. */

	/*  N       (input) INTEGER */
	/*          The number of columns of the matrix A.  N >= 0. */

	/*  A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          The m by n matrix A.  If UPLO = 'U', only the upper triangle */
	/*          or trapezoid is accessed; if UPLO = 'L', only the lower */
	/*          triangle or trapezoid is accessed. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(1,M). */

	/*  B       (output) DOUBLE PRECISION array, dimension (LDB,N) */
	/*          On exit, B = A in the locations specified by UPLO. */

	/*  LDB     (input) INTEGER */
	/*          The leading dimension of the array B.  LDB >= max(1,M). */

	/*  ===================================================================== */

	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	b_dim1 = *ldb;
	b_offset = 1 + b_dim1;
	b -= b_offset;

	/* Function Body */
	if (lsame_ (uplo, "U", (ftnlen) 1, (ftnlen) 1))
	{
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			i__2 = min (j, *m);
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				b[i__ + j * b_dim1] = a[i__ + j * a_dim1];
				/* L10: */
			}
			/* L20: */
		}
	}
	else if (lsame_ (uplo, "L", (ftnlen) 1, (ftnlen) 1))
	{
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			i__2 = *m;
			for (i__ = j; i__ <= i__2; ++i__)
			{
				b[i__ + j * b_dim1] = a[i__ + j * a_dim1];
				/* L30: */
			}
			/* L40: */
		}
	}
	else
	{
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				b[i__ + j * b_dim1] = a[i__ + j * a_dim1];
				/* L50: */
			}
			/* L60: */
		}
	}
	return 0;

	/*     End of DLACPY */

}	/* dlacpy_ */

static int
dladiv_ (doublereal * a, doublereal * b, doublereal * c__, doublereal * d__, doublereal * p, doublereal * q)
{
	static doublereal e, f;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     October 31, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLADIV performs complex division in  real arithmetic */

	/*                        a + i*b */
	/*             p + i*q = --------- */
	/*                        c + i*d */

	/*  The algorithm is due to Robert L. Smith and can be found */
	/*  in D. Knuth, The art of Computer Programming, Vol.2, p.195 */

	/*  Arguments */
	/*  ========= */

	/*  A       (input) DOUBLE PRECISION */
	/*  B       (input) DOUBLE PRECISION */
	/*  C       (input) DOUBLE PRECISION */
	/*  D       (input) DOUBLE PRECISION */
	/*          The scalars a, b, c, and d in the above expression. */

	/*  P       (output) DOUBLE PRECISION */
	/*  Q       (output) DOUBLE PRECISION */
	/*          The scalars p and q in the above expression. */

	/*  ===================================================================== */

	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	if (abs (*d__) < abs (*c__))
	{
		e = *d__ / *c__;
		f = *c__ + *d__ * e;
		*p = (*a + *b * e) / f;
		*q = (*b - *a * e) / f;
	}
	else
	{
		e = *c__ / *d__;
		f = *d__ + *c__ * e;
		*p = (*b + *a * e) / f;
		*q = (-(*a) + *b * e) / f;
	}

	return 0;

	/*     End of DLADIV */

}	/* dladiv_ */

static int
dlaexc_ (logical * wantq, integer * n, doublereal * t,
			integer * ldt, doublereal * q, integer * ldq, integer * j1, integer * n1, integer * n2, doublereal * work, integer * info)
{
	/* System generated locals */
	integer q_dim1, q_offset, t_dim1, t_offset, i__1;
	doublereal d__1, d__2, d__3;

	/* Local variables */
	static doublereal d__[16] /* was [4][4] */ ;
	static integer k;
	static doublereal u[3], x[4] /* was [2][2] */ ;
	static integer j2, j3, j4;
	static doublereal u1[3], u2[3];
	static integer nd;
	static doublereal cs, t11, t22, t33, sn, wi1, wi2, wr1, wr2, eps, tau, tau1, tau2;
	static integer ierr;
	static doublereal temp;
	static doublereal scale, dnorm, xnorm;
	static doublereal thresh, smlnum;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     February 29, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLAEXC swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in */
	/*  an upper quasi-triangular matrix T by an orthogonal similarity */
	/*  transformation. */

	/*  T must be in Schur canonical form, that is, block upper triangular */
	/*  with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block */
	/*  has its diagonal elemnts equal and its off-diagonal elements of */
	/*  opposite sign. */

	/*  Arguments */
	/*  ========= */

	/*  WANTQ   (input) LOGICAL */
	/*          = .TRUE. : accumulate the transformation in the matrix Q; */
	/*          = .FALSE.: do not accumulate the transformation. */

	/*  N       (input) INTEGER */
	/*          The order of the matrix T. N >= 0. */

	/*  T       (input/output) DOUBLE PRECISION array, dimension (LDT,N) */
	/*          On entry, the upper quasi-triangular matrix T, in Schur */
	/*          canonical form. */
	/*          On exit, the updated matrix T, again in Schur canonical form. */

	/*  LDT     (input)  INTEGER */
	/*          The leading dimension of the array T. LDT >= max(1,N). */

	/*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
	/*          On entry, if WANTQ is .TRUE., the orthogonal matrix Q. */
	/*          On exit, if WANTQ is .TRUE., the updated matrix Q. */
	/*          If WANTQ is .FALSE., Q is not referenced. */

	/*  LDQ     (input) INTEGER */
	/*          The leading dimension of the array Q. */
	/*          LDQ >= 1; and if WANTQ is .TRUE., LDQ >= N. */

	/*  J1      (input) INTEGER */
	/*          The index of the first row of the first block T11. */

	/*  N1      (input) INTEGER */
	/*          The order of the first block T11. N1 = 0, 1 or 2. */

	/*  N2      (input) INTEGER */
	/*          The order of the second block T22. N2 = 0, 1 or 2. */

	/*  WORK    (workspace) DOUBLE PRECISION array, dimension (N) */

	/*  INFO    (output) INTEGER */
	/*          = 0: successful exit */
	/*          = 1: the transformed matrix T would be too far from Schur */
	/*               form; the blocks are not swapped and T and Q are */
	/*               unchanged. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Local Arrays .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/* Parameter adjustments */
	t_dim1 = *ldt;
	t_offset = 1 + t_dim1;
	t -= t_offset;
	q_dim1 = *ldq;
	q_offset = 1 + q_dim1;
	q -= q_offset;
	--work;

	/* Function Body */
	*info = 0;

	/*     Quick return if possible */

	if (*n == 0 || *n1 == 0 || *n2 == 0)
	{
		return 0;
	}
	if (*j1 + *n1 > *n)
	{
		return 0;
	}

	j2 = *j1 + 1;
	j3 = *j1 + 2;
	j4 = *j1 + 3;

	if (*n1 == 1 && *n2 == 1)
	{

		/*        Swap two 1-by-1 blocks. */

		t11 = t[*j1 + *j1 * t_dim1];
		t22 = t[j2 + j2 * t_dim1];

		/*        Determine the transformation to perform the interchange. */

		d__1 = t22 - t11;
		dlartg_ (&t[*j1 + j2 * t_dim1], &d__1, &cs, &sn, &temp);

		/*        Apply transformation to the matrix T. */

		if (j3 <= *n)
		{
			i__1 = *n - *j1 - 1;
			drot_ (&i__1, &t[*j1 + j3 * t_dim1], ldt, &t[j2 + j3 * t_dim1], ldt, &cs, &sn);
		}
		i__1 = *j1 - 1;
		drot_ (&i__1, &t[*j1 * t_dim1 + 1], &c__1, &t[j2 * t_dim1 + 1], &c__1, &cs, &sn);

		t[*j1 + *j1 * t_dim1] = t22;
		t[j2 + j2 * t_dim1] = t11;

		if (*wantq)
		{

			/*           Accumulate transformation in the matrix Q. */

			drot_ (n, &q[*j1 * q_dim1 + 1], &c__1, &q[j2 * q_dim1 + 1], &c__1, &cs, &sn);
		}

	}
	else
	{

		/*        Swapping involves at least one 2-by-2 block. */

		/*        Copy the diagonal block of order N1+N2 to the local array D */
		/*        and compute its norm. */

		nd = *n1 + *n2;
		dlacpy_ ("Full", &nd, &nd, &t[*j1 + *j1 * t_dim1], ldt, d__, &c__4, (ftnlen) 4);
		dnorm = dlange_ ("Max", &nd, &nd, d__, &c__4, &work[1], (ftnlen) 3);

		/*        Compute machine-dependent threshold for test for accepting */
		/*        swap. */

		eps = dlamch_ ("P", (ftnlen) 1);
		smlnum = dlamch_ ("S", (ftnlen) 1) / eps;
		/* Computing MAX */
		d__1 = eps * 10. * dnorm;
		thresh = max (d__1, smlnum);

		/*        Solve T11*X - X*T22 = scale*T12 for X. */

		dlasy2_ (&c_false, &c_false, &c_n1, n1, n2, d__, &c__4, &d__[*n1 + 1 +
																						 (*n1 + 1 << 2) - 5], &c__4, &d__[(*n1 + 1 << 2) - 4], &c__4, &scale, x, &c__2, &xnorm,
					&ierr);

		/*        Swap the adjacent diagonal blocks. */

		k = *n1 + *n1 + *n2 - 3;
		switch (k)
		{
		case 1:
			goto L10;
		case 2:
			goto L20;
		case 3:
			goto L30;
		}

	 L10:

		/*        N1 = 1, N2 = 2: generate elementary reflector H so that: */

		/*        ( scale, X11, X12 ) H = ( 0, 0, * ) */

		u[0] = scale;
		u[1] = x[0];
		u[2] = x[2];
		dlarfg_ (&c__3, &u[2], u, &c__1, &tau);
		u[2] = 1.;
		t11 = t[*j1 + *j1 * t_dim1];

		/*        Perform swap provisionally on diagonal block in D. */

		dlarfx_ ("L", &c__3, &c__3, u, &tau, d__, &c__4, &work[1], (ftnlen) 1);
		dlarfx_ ("R", &c__3, &c__3, u, &tau, d__, &c__4, &work[1], (ftnlen) 1);

		/*        Test whether to reject swap. */

		/* Computing MAX */
		d__2 = abs (d__[2]), d__3 = abs (d__[6]), d__2 = max (d__2, d__3), d__3 = (d__1 = d__[10] - t11, abs (d__1));
		if (max (d__2, d__3) > thresh)
		{
			goto L50;
		}

		/*        Accept swap: apply transformation to the entire matrix T. */

		i__1 = *n - *j1 + 1;
		dlarfx_ ("L", &c__3, &i__1, u, &tau, &t[*j1 + *j1 * t_dim1], ldt, &work[1], (ftnlen) 1);
		dlarfx_ ("R", &j2, &c__3, u, &tau, &t[*j1 * t_dim1 + 1], ldt, &work[1], (ftnlen) 1);

		t[j3 + *j1 * t_dim1] = 0.;
		t[j3 + j2 * t_dim1] = 0.;
		t[j3 + j3 * t_dim1] = t11;

		if (*wantq)
		{

			/*           Accumulate transformation in the matrix Q. */

			dlarfx_ ("R", n, &c__3, u, &tau, &q[*j1 * q_dim1 + 1], ldq, &work[1], (ftnlen) 1);
		}
		goto L40;

	 L20:

		/*        N1 = 2, N2 = 1: generate elementary reflector H so that: */

		/*        H (  -X11 ) = ( * ) */
		/*          (  -X21 ) = ( 0 ) */
		/*          ( scale ) = ( 0 ) */

		u[0] = -x[0];
		u[1] = -x[1];
		u[2] = scale;
		dlarfg_ (&c__3, u, &u[1], &c__1, &tau);
		u[0] = 1.;
		t33 = t[j3 + j3 * t_dim1];

		/*        Perform swap provisionally on diagonal block in D. */

		dlarfx_ ("L", &c__3, &c__3, u, &tau, d__, &c__4, &work[1], (ftnlen) 1);
		dlarfx_ ("R", &c__3, &c__3, u, &tau, d__, &c__4, &work[1], (ftnlen) 1);

		/*        Test whether to reject swap. */

		/* Computing MAX */
		d__2 = abs (d__[1]), d__3 = abs (d__[2]), d__2 = max (d__2, d__3), d__3 = (d__1 = d__[0] - t33, abs (d__1));
		if (max (d__2, d__3) > thresh)
		{
			goto L50;
		}

		/*        Accept swap: apply transformation to the entire matrix T. */

		dlarfx_ ("R", &j3, &c__3, u, &tau, &t[*j1 * t_dim1 + 1], ldt, &work[1], (ftnlen) 1);
		i__1 = *n - *j1;
		dlarfx_ ("L", &c__3, &i__1, u, &tau, &t[*j1 + j2 * t_dim1], ldt, &work[1], (ftnlen) 1);

		t[*j1 + *j1 * t_dim1] = t33;
		t[j2 + *j1 * t_dim1] = 0.;
		t[j3 + *j1 * t_dim1] = 0.;

		if (*wantq)
		{

			/*           Accumulate transformation in the matrix Q. */

			dlarfx_ ("R", n, &c__3, u, &tau, &q[*j1 * q_dim1 + 1], ldq, &work[1], (ftnlen) 1);
		}
		goto L40;

	 L30:

		/*        N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so */
		/*        that: */

		/*        H(2) H(1) (  -X11  -X12 ) = (  *  * ) */
		/*                  (  -X21  -X22 )   (  0  * ) */
		/*                  ( scale    0  )   (  0  0 ) */
		/*                  (    0  scale )   (  0  0 ) */

		u1[0] = -x[0];
		u1[1] = -x[1];
		u1[2] = scale;
		dlarfg_ (&c__3, u1, &u1[1], &c__1, &tau1);
		u1[0] = 1.;

		temp = -tau1 * (x[2] + u1[1] * x[3]);
		u2[0] = -temp * u1[1] - x[3];
		u2[1] = -temp * u1[2];
		u2[2] = scale;
		dlarfg_ (&c__3, u2, &u2[1], &c__1, &tau2);
		u2[0] = 1.;

		/*        Perform swap provisionally on diagonal block in D. */

		dlarfx_ ("L", &c__3, &c__4, u1, &tau1, d__, &c__4, &work[1], (ftnlen) 1);
		dlarfx_ ("R", &c__4, &c__3, u1, &tau1, d__, &c__4, &work[1], (ftnlen) 1);
		dlarfx_ ("L", &c__3, &c__4, u2, &tau2, &d__[1], &c__4, &work[1], (ftnlen) 1);
		dlarfx_ ("R", &c__4, &c__3, u2, &tau2, &d__[4], &c__4, &work[1], (ftnlen) 1);

		/*        Test whether to reject swap. */

		/* Computing MAX */
		d__1 = abs (d__[2]), d__2 = abs (d__[6]), d__1 = max (d__1, d__2), d__2 = abs (d__[3]), d__1 = max (d__1, d__2), d__2 = abs (d__[7]);
		if (max (d__1, d__2) > thresh)
		{
			goto L50;
		}

		/*        Accept swap: apply transformation to the entire matrix T. */

		i__1 = *n - *j1 + 1;
		dlarfx_ ("L", &c__3, &i__1, u1, &tau1, &t[*j1 + *j1 * t_dim1], ldt, &work[1], (ftnlen) 1);
		dlarfx_ ("R", &j4, &c__3, u1, &tau1, &t[*j1 * t_dim1 + 1], ldt, &work[1], (ftnlen) 1);
		i__1 = *n - *j1 + 1;
		dlarfx_ ("L", &c__3, &i__1, u2, &tau2, &t[j2 + *j1 * t_dim1], ldt, &work[1], (ftnlen) 1);
		dlarfx_ ("R", &j4, &c__3, u2, &tau2, &t[j2 * t_dim1 + 1], ldt, &work[1], (ftnlen) 1);

		t[j3 + *j1 * t_dim1] = 0.;
		t[j3 + j2 * t_dim1] = 0.;
		t[j4 + *j1 * t_dim1] = 0.;
		t[j4 + j2 * t_dim1] = 0.;

		if (*wantq)
		{

			/*           Accumulate transformation in the matrix Q. */

			dlarfx_ ("R", n, &c__3, u1, &tau1, &q[*j1 * q_dim1 + 1], ldq, &work[1], (ftnlen) 1);
			dlarfx_ ("R", n, &c__3, u2, &tau2, &q[j2 * q_dim1 + 1], ldq, &work[1], (ftnlen) 1);
		}

	 L40:

		if (*n2 == 2)
		{

			/*           Standardize new 2-by-2 block T11 */

			dlanv2_ (&t[*j1 + *j1 * t_dim1], &t[*j1 + j2 * t_dim1], &t[j2 + *j1 * t_dim1], &t[j2 + j2 * t_dim1], &wr1, &wi1, &wr2, &wi2, &cs, &sn);
			i__1 = *n - *j1 - 1;
			drot_ (&i__1, &t[*j1 + (*j1 + 2) * t_dim1], ldt, &t[j2 + (*j1 + 2) * t_dim1], ldt, &cs, &sn);
			i__1 = *j1 - 1;
			drot_ (&i__1, &t[*j1 * t_dim1 + 1], &c__1, &t[j2 * t_dim1 + 1], &c__1, &cs, &sn);
			if (*wantq)
			{
				drot_ (n, &q[*j1 * q_dim1 + 1], &c__1, &q[j2 * q_dim1 + 1], &c__1, &cs, &sn);
			}
		}

		if (*n1 == 2)
		{

			/*           Standardize new 2-by-2 block T22 */

			j3 = *j1 + *n2;
			j4 = j3 + 1;
			dlanv2_ (&t[j3 + j3 * t_dim1], &t[j3 + j4 * t_dim1], &t[j4 + j3 * t_dim1], &t[j4 + j4 * t_dim1], &wr1, &wi1, &wr2, &wi2, &cs, &sn);
			if (j3 + 2 <= *n)
			{
				i__1 = *n - j3 - 1;
				drot_ (&i__1, &t[j3 + (j3 + 2) * t_dim1], ldt, &t[j4 + (j3 + 2) * t_dim1], ldt, &cs, &sn);
			}
			i__1 = j3 - 1;
			drot_ (&i__1, &t[j3 * t_dim1 + 1], &c__1, &t[j4 * t_dim1 + 1], &c__1, &cs, &sn);
			if (*wantq)
			{
				drot_ (n, &q[j3 * q_dim1 + 1], &c__1, &q[j4 * q_dim1 + 1], &c__1, &cs, &sn);
			}
		}

	}
	return 0;

	/*     Exit with INFO = 1 if swap was rejected. */

 L50:
	*info = 1;
	return 0;

	/*     End of DLAEXC */

}	/* dlaexc_ */

static int
dlahqr_ (logical * wantt, logical * wantz, integer * n,
			integer * ilo, integer * ihi, doublereal * h__, integer * ldh, doublereal
			* wr, doublereal * wi, integer * iloz, integer * ihiz, doublereal * z__, integer * ldz, integer * info)
{
	/* System generated locals */
	integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4;
	doublereal d__1, d__2;

	/* Builtin functions */
	double sqrt (doublereal), d_sign (doublereal *, doublereal *);

	/* Local variables */
	static integer i__, j, k, l, m;
	static doublereal s, v[3];
	static integer i1, i2;
	static doublereal t1, t2, t3, v1, v2, v3, h00, h10, h11, h12, h21, h22, h33, h44;
	static integer nh;
	static doublereal cs;
	static integer nr;
	static doublereal sn;
	static integer nz;
	static doublereal ave, h33s, h44s;
	static integer itn, its;
	static doublereal ulp, sum, tst1, h43h34, disc, unfl, ovfl;
	static doublereal work[1];
	static doublereal smlnum;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     June 30, 1999 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLAHQR is an auxiliary routine called by DHSEQR to update the */
	/*  eigenvalues and Schur decomposition already computed by DHSEQR, by */
	/*  dealing with the Hessenberg submatrix in rows and columns ILO to IHI. */

	/*  Arguments */
	/*  ========= */

	/*  WANTT   (input) LOGICAL */
	/*          = .TRUE. : the full Schur form T is required; */
	/*          = .FALSE.: only eigenvalues are required. */

	/*  WANTZ   (input) LOGICAL */
	/*          = .TRUE. : the matrix of Schur vectors Z is required; */
	/*          = .FALSE.: Schur vectors are not required. */

	/*  N       (input) INTEGER */
	/*          The order of the matrix H.  N >= 0. */

	/*  ILO     (input) INTEGER */
	/*  IHI     (input) INTEGER */
	/*          It is assumed that H is already upper quasi-triangular in */
	/*          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless */
	/*          ILO = 1). DLAHQR works primarily with the Hessenberg */
	/*          submatrix in rows and columns ILO to IHI, but applies */
	/*          transformations to all of H if WANTT is .TRUE.. */
	/*          1 <= ILO <= max(1,IHI); IHI <= N. */

	/*  H       (input/output) DOUBLE PRECISION array, dimension (LDH,N) */
	/*          On entry, the upper Hessenberg matrix H. */
	/*          On exit, if WANTT is .TRUE., H is upper quasi-triangular in */
	/*          rows and columns ILO:IHI, with any 2-by-2 diagonal blocks in */
	/*          standard form. If WANTT is .FALSE., the contents of H are */
	/*          unspecified on exit. */

	/*  LDH     (input) INTEGER */
	/*          The leading dimension of the array H. LDH >= max(1,N). */

	/*  WR      (output) DOUBLE PRECISION array, dimension (N) */
	/*  WI      (output) DOUBLE PRECISION array, dimension (N) */
	/*          The real and imaginary parts, respectively, of the computed */
	/*          eigenvalues ILO to IHI are stored in the corresponding */
	/*          elements of WR and WI. If two eigenvalues are computed as a */
	/*          complex conjugate pair, they are stored in consecutive */
	/*          elements of WR and WI, say the i-th and (i+1)th, with */
	/*          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the */
	/*          eigenvalues are stored in the same order as on the diagonal */
	/*          of the Schur form returned in H, with WR(i) = H(i,i), and, if */
	/*          H(i:i+1,i:i+1) is a 2-by-2 diagonal block, */
	/*          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i). */

	/*  ILOZ    (input) INTEGER */
	/*  IHIZ    (input) INTEGER */
	/*          Specify the rows of Z to which transformations must be */
	/*          applied if WANTZ is .TRUE.. */
	/*          1 <= ILOZ <= ILO; IHI <= IHIZ <= N. */

	/*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N) */
	/*          If WANTZ is .TRUE., on entry Z must contain the current */
	/*          matrix Z of transformations accumulated by DHSEQR, and on */
	/*          exit Z has been updated; transformations are applied only to */
	/*          the submatrix Z(ILOZ:IHIZ,ILO:IHI). */
	/*          If WANTZ is .FALSE., Z is not referenced. */

	/*  LDZ     (input) INTEGER */
	/*          The leading dimension of the array Z. LDZ >= max(1,N). */

	/*  INFO    (output) INTEGER */
	/*          = 0: successful exit */
	/*          > 0: DLAHQR failed to compute all the eigenvalues ILO to IHI */
	/*               in a total of 30*(IHI-ILO+1) iterations; if INFO = i, */
	/*               elements i+1:ihi of WR and WI contain those eigenvalues */
	/*               which have been successfully computed. */

	/*  Further Details */
	/*  =============== */

	/*  2-96 Based on modifications by */
	/*     David Day, Sandia National Laboratory, USA */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Local Arrays .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/* Parameter adjustments */
	h_dim1 = *ldh;
	h_offset = 1 + h_dim1;
	h__ -= h_offset;
	--wr;
	--wi;
	z_dim1 = *ldz;
	z_offset = 1 + z_dim1;
	z__ -= z_offset;

	/* Function Body */
	*info = 0;

	/*     Quick return if possible */

	if (*n == 0)
	{
		return 0;
	}
	if (*ilo == *ihi)
	{
		wr[*ilo] = h__[*ilo + *ilo * h_dim1];
		wi[*ilo] = 0.;
		return 0;
	}

	nh = *ihi - *ilo + 1;
	nz = *ihiz - *iloz + 1;

	/*     Set machine-dependent constants for the stopping criterion. */
	/*     If norm(H) <= sqrt(OVFL), overflow should not occur. */

	unfl = dlamch_ ("Safe minimum", (ftnlen) 12);
	ovfl = 1. / unfl;
	dlabad_ (&unfl, &ovfl);
	ulp = dlamch_ ("Precision", (ftnlen) 9);
	smlnum = unfl * (nh / ulp);

	/*     I1 and I2 are the indices of the first row and last column of H */
	/*     to which transformations must be applied. If eigenvalues only are */
	/*     being computed, I1 and I2 are set inside the main loop. */

	if (*wantt)
	{
		i1 = 1;
		i2 = *n;
	}

	/*     ITN is the total number of QR iterations allowed. */

	itn = nh * 30;

	/*     The main loop begins here. I is the loop index and decreases from */
	/*     IHI to ILO in steps of 1 or 2. Each iteration of the loop works */
	/*     with the active submatrix in rows and columns L to I. */
	/*     Eigenvalues I+1 to IHI have already converged. Either L = ILO or */
	/*     H(L,L-1) is negligible so that the matrix splits. */

	i__ = *ihi;
 L10:
	l = *ilo;
	if (i__ < *ilo)
	{
		goto L150;
	}

	/*     Perform QR iterations on rows and columns ILO to I until a */
	/*     submatrix of order 1 or 2 splits off at the bottom because a */
	/*     subdiagonal element has become negligible. */

	i__1 = itn;
	for (its = 0; its <= i__1; ++its)
	{

		/*        Look for a single small subdiagonal element. */

		i__2 = l + 1;
		for (k = i__; k >= i__2; --k)
		{
			tst1 = (d__1 = h__[k - 1 + (k - 1) * h_dim1], abs (d__1)) + (d__2 = h__[k + k * h_dim1], abs (d__2));
			if (tst1 == 0.)
			{
				i__3 = i__ - l + 1;
				tst1 = dlanhs_ ("1", &i__3, &h__[l + l * h_dim1], ldh, work, (ftnlen) 1);
			}
			/* Computing MAX */
			d__2 = ulp * tst1;
			if ((d__1 = h__[k + (k - 1) * h_dim1], abs (d__1)) <= max (d__2, smlnum))
			{
				goto L30;
			}
			/* L20: */
		}
	 L30:
		l = k;
		if (l > *ilo)
		{

			/*           H(L,L-1) is negligible */

			h__[l + (l - 1) * h_dim1] = 0.;
		}

		/*        Exit from loop if a submatrix of order 1 or 2 has split off. */

		if (l >= i__ - 1)
		{
			goto L140;
		}

		/*        Now the active submatrix is in rows and columns L to I. If */
		/*        eigenvalues only are being computed, only the active submatrix */
		/*        need be transformed. */

		if (!(*wantt))
		{
			i1 = l;
			i2 = i__;
		}

		if (its == 10 || its == 20)
		{

			/*           Exceptional shift. */

			s = (d__1 = h__[i__ + (i__ - 1) * h_dim1], abs (d__1)) + (d__2 = h__[i__ - 1 + (i__ - 2) * h_dim1], abs (d__2));
			h44 = s * .75 + h__[i__ + i__ * h_dim1];
			h33 = h44;
			h43h34 = s * -.4375 * s;
		}
		else
		{

			/*           Prepare to use Francis' double shift */
			/*           (i.e. 2nd degree generalized Rayleigh quotient) */

			h44 = h__[i__ + i__ * h_dim1];
			h33 = h__[i__ - 1 + (i__ - 1) * h_dim1];
			h43h34 = h__[i__ + (i__ - 1) * h_dim1] * h__[i__ - 1 + i__ * h_dim1];
			s = h__[i__ - 1 + (i__ - 2) * h_dim1] * h__[i__ - 1 + (i__ - 2) * h_dim1];
			disc = (h33 - h44) * .5;
			disc = disc * disc + h43h34;
			if (disc > 0.)
			{

				/*              Real roots: use Wilkinson's shift twice */

				disc = sqrt (disc);
				ave = (h33 + h44) * .5;
				if (abs (h33) - abs (h44) > 0.)
				{
					h33 = h33 * h44 - h43h34;
					h44 = h33 / (d_sign (&disc, &ave) + ave);
				}
				else
				{
					h44 = d_sign (&disc, &ave) + ave;
				}
				h33 = h44;
				h43h34 = 0.;
			}
		}

		/*        Look for two consecutive small subdiagonal elements. */

		i__2 = l;
		for (m = i__ - 2; m >= i__2; --m)
		{
			/*           Determine the effect of starting the double-shift QR */
			/*           iteration at row M, and see if this would make H(M,M-1) */
			/*           negligible. */

			h11 = h__[m + m * h_dim1];
			h22 = h__[m + 1 + (m + 1) * h_dim1];
			h21 = h__[m + 1 + m * h_dim1];
			h12 = h__[m + (m + 1) * h_dim1];
			h44s = h44 - h11;
			h33s = h33 - h11;
			v1 = (h33s * h44s - h43h34) / h21 + h12;
			v2 = h22 - h11 - h33s - h44s;
			v3 = h__[m + 2 + (m + 1) * h_dim1];
			s = abs (v1) + abs (v2) + abs (v3);
			v1 /= s;
			v2 /= s;
			v3 /= s;
			v[0] = v1;
			v[1] = v2;
			v[2] = v3;
			if (m == l)
			{
				goto L50;
			}
			h00 = h__[m - 1 + (m - 1) * h_dim1];
			h10 = h__[m + (m - 1) * h_dim1];
			tst1 = abs (v1) * (abs (h00) + abs (h11) + abs (h22));
			if (abs (h10) * (abs (v2) + abs (v3)) <= ulp * tst1)
			{
				goto L50;
			}
			/* L40: */
		}
	 L50:

		/*        Double-shift QR step */

		i__2 = i__ - 1;
		for (k = m; k <= i__2; ++k)
		{

			/*           The first iteration of this loop determines a reflection G */
			/*           from the vector V and applies it from left and right to H, */
			/*           thus creating a nonzero bulge below the subdiagonal. */

			/*           Each subsequent iteration determines a reflection G to */
			/*           restore the Hessenberg form in the (K-1)th column, and thus */
			/*           chases the bulge one step toward the bottom of the active */
			/*           submatrix. NR is the order of G. */

			/* Computing MIN */
			i__3 = 3, i__4 = i__ - k + 1;
			nr = min (i__3, i__4);
			if (k > m)
			{
				dcopy_ (&nr, &h__[k + (k - 1) * h_dim1], &c__1, v, &c__1);
			}
			dlarfg_ (&nr, v, &v[1], &c__1, &t1);
			if (k > m)
			{
				h__[k + (k - 1) * h_dim1] = v[0];
				h__[k + 1 + (k - 1) * h_dim1] = 0.;
				if (k < i__ - 1)
				{
					h__[k + 2 + (k - 1) * h_dim1] = 0.;
				}
			}
			else if (m > l)
			{
				h__[k + (k - 1) * h_dim1] = -h__[k + (k - 1) * h_dim1];
			}
			v2 = v[1];
			t2 = t1 * v2;
			if (nr == 3)
			{
				v3 = v[2];
				t3 = t1 * v3;

				/*              Apply G from the left to transform the rows of the matrix */
				/*              in columns K to I2. */

				i__3 = i2;
				for (j = k; j <= i__3; ++j)
				{
					sum = h__[k + j * h_dim1] + v2 * h__[k + 1 + j * h_dim1] + v3 * h__[k + 2 + j * h_dim1];
					h__[k + j * h_dim1] -= sum * t1;
					h__[k + 1 + j * h_dim1] -= sum * t2;
					h__[k + 2 + j * h_dim1] -= sum * t3;
					/* L60: */
				}

				/*              Apply G from the right to transform the columns of the */
				/*              matrix in rows I1 to min(K+3,I). */

				/* Computing MIN */
				i__4 = k + 3;
				i__3 = min (i__4, i__);
				for (j = i1; j <= i__3; ++j)
				{
					sum = h__[j + k * h_dim1] + v2 * h__[j + (k + 1) * h_dim1] + v3 * h__[j + (k + 2) * h_dim1];
					h__[j + k * h_dim1] -= sum * t1;
					h__[j + (k + 1) * h_dim1] -= sum * t2;
					h__[j + (k + 2) * h_dim1] -= sum * t3;
					/* L70: */
				}

				if (*wantz)
				{

					/*                 Accumulate transformations in the matrix Z */

					i__3 = *ihiz;
					for (j = *iloz; j <= i__3; ++j)
					{
						sum = z__[j + k * z_dim1] + v2 * z__[j + (k + 1) * z_dim1] + v3 * z__[j + (k + 2) * z_dim1];
						z__[j + k * z_dim1] -= sum * t1;
						z__[j + (k + 1) * z_dim1] -= sum * t2;
						z__[j + (k + 2) * z_dim1] -= sum * t3;
						/* L80: */
					}
				}
			}
			else if (nr == 2)
			{

				/*              Apply G from the left to transform the rows of the matrix */
				/*              in columns K to I2. */

				i__3 = i2;
				for (j = k; j <= i__3; ++j)
				{
					sum = h__[k + j * h_dim1] + v2 * h__[k + 1 + j * h_dim1];
					h__[k + j * h_dim1] -= sum * t1;
					h__[k + 1 + j * h_dim1] -= sum * t2;
					/* L90: */
				}

				/*              Apply G from the right to transform the columns of the */
				/*              matrix in rows I1 to min(K+3,I). */

				i__3 = i__;
				for (j = i1; j <= i__3; ++j)
				{
					sum = h__[j + k * h_dim1] + v2 * h__[j + (k + 1) * h_dim1];
					h__[j + k * h_dim1] -= sum * t1;
					h__[j + (k + 1) * h_dim1] -= sum * t2;
					/* L100: */
				}

				if (*wantz)
				{

					/*                 Accumulate transformations in the matrix Z */

					i__3 = *ihiz;
					for (j = *iloz; j <= i__3; ++j)
					{
						sum = z__[j + k * z_dim1] + v2 * z__[j + (k + 1) * z_dim1];
						z__[j + k * z_dim1] -= sum * t1;
						z__[j + (k + 1) * z_dim1] -= sum * t2;
						/* L110: */
					}
				}
			}
			/* L120: */
		}

		/* L130: */
	}

	/*     Failure to converge in remaining number of iterations */

	*info = i__;
	return 0;

 L140:

	if (l == i__)
	{

		/*        H(I,I-1) is negligible: one eigenvalue has converged. */

		wr[i__] = h__[i__ + i__ * h_dim1];
		wi[i__] = 0.;
	}
	else if (l == i__ - 1)
	{

		/*        H(I-1,I-2) is negligible: a pair of eigenvalues have converged. */

		/*        Transform the 2-by-2 submatrix to standard Schur form, */
		/*        and compute and store the eigenvalues. */

		dlanv2_ (&h__[i__ - 1 + (i__ - 1) * h_dim1], &h__[i__ - 1 + i__ *
																		  h_dim1], &h__[i__ + (i__ - 1) * h_dim1], &h__[i__ + i__ *
																																		h_dim1], &wr[i__ - 1], &wi[i__ - 1], &wr[i__], &wi[i__],
					&cs, &sn);

		if (*wantt)
		{

			/*           Apply the transformation to the rest of H. */

			if (i2 > i__)
			{
				i__1 = i2 - i__;
				drot_ (&i__1, &h__[i__ - 1 + (i__ + 1) * h_dim1], ldh, &h__[i__ + (i__ + 1) * h_dim1], ldh, &cs, &sn);
			}
			i__1 = i__ - i1 - 1;
			drot_ (&i__1, &h__[i1 + (i__ - 1) * h_dim1], &c__1, &h__[i1 + i__ * h_dim1], &c__1, &cs, &sn);
		}
		if (*wantz)
		{

			/*           Apply the transformation to Z. */

			drot_ (&nz, &z__[*iloz + (i__ - 1) * z_dim1], &c__1, &z__[*iloz + i__ * z_dim1], &c__1, &cs, &sn);
		}
	}

	/*     Decrement number of remaining iterations, and return to start of */
	/*     the main loop with new value of I. */

	itn -= its;
	i__ = l - 1;
	goto L10;

 L150:
	return 0;

	/*     End of DLAHQR */

}	/* dlahqr_ */

static int
dlahrd_ (integer * n, integer * k, integer * nb, doublereal * a, integer * lda, doublereal * tau, doublereal * t, integer * ldt, doublereal * y, integer * ldy)
{
	/* System generated locals */
	integer a_dim1, a_offset, t_dim1, t_offset, y_dim1, y_offset, i__1, i__2, i__3;
	doublereal d__1;

	/* Local variables */
	static integer i__;
	static doublereal ei;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     June 30, 1999 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLAHRD reduces the first NB columns of a real general n-by-(n-k+1) */
	/*  matrix A so that elements below the k-th subdiagonal are zero. The */
	/*  reduction is performed by an orthogonal similarity transformation */
	/*  Q' * A * Q. The routine returns the matrices V and T which determine */
	/*  Q as a block reflector I - V*T*V', and also the matrix Y = A * V * T. */

	/*  This is an auxiliary routine called by DGEHRD. */

	/*  Arguments */
	/*  ========= */

	/*  N       (input) INTEGER */
	/*          The order of the matrix A. */

	/*  K       (input) INTEGER */
	/*          The offset for the reduction. Elements below the k-th */
	/*          subdiagonal in the first NB columns are reduced to zero. */

	/*  NB      (input) INTEGER */
	/*          The number of columns to be reduced. */

	/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N-K+1) */
	/*          On entry, the n-by-(n-k+1) general matrix A. */
	/*          On exit, the elements on and above the k-th subdiagonal in */
	/*          the first NB columns are overwritten with the corresponding */
	/*          elements of the reduced matrix; the elements below the k-th */
	/*          subdiagonal, with the array TAU, represent the matrix Q as a */
	/*          product of elementary reflectors. The other columns of A are */
	/*          unchanged. See Further Details. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(1,N). */

	/*  TAU     (output) DOUBLE PRECISION array, dimension (NB) */
	/*          The scalar factors of the elementary reflectors. See Further */
	/*          Details. */

	/*  T       (output) DOUBLE PRECISION array, dimension (LDT,NB) */
	/*          The upper triangular matrix T. */

	/*  LDT     (input) INTEGER */
	/*          The leading dimension of the array T.  LDT >= NB. */

	/*  Y       (output) DOUBLE PRECISION array, dimension (LDY,NB) */
	/*          The n-by-nb matrix Y. */

	/*  LDY     (input) INTEGER */
	/*          The leading dimension of the array Y. LDY >= N. */

	/*  Further Details */
	/*  =============== */

	/*  The matrix Q is represented as a product of nb elementary reflectors */

	/*     Q = H(1) H(2) . . . H(nb). */

	/*  Each H(i) has the form */

	/*     H(i) = I - tau * v * v' */

	/*  where tau is a real scalar, and v is a real vector with */
	/*  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in */
	/*  A(i+k+1:n,i), and tau in TAU(i). */

	/*  The elements of the vectors v together form the (n-k+1)-by-nb matrix */
	/*  V which is needed, with T and Y, to apply the transformation to the */
	/*  unreduced part of the matrix, using an update of the form: */
	/*  A := (I - V*T*V') * (A - Y*V'). */

	/*  The contents of A on exit are illustrated by the following example */
	/*  with n = 7, k = 3 and nb = 2: */

	/*     ( a   h   a   a   a ) */
	/*     ( a   h   a   a   a ) */
	/*     ( a   h   a   a   a ) */
	/*     ( h   h   a   a   a ) */
	/*     ( v1  h   a   a   a ) */
	/*     ( v1  v2  a   a   a ) */
	/*     ( v1  v2  a   a   a ) */

	/*  where a denotes an element of the original matrix A, h denotes a */
	/*  modified element of the upper Hessenberg matrix H, and vi denotes an */
	/*  element of the vector defining H(i). */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Quick return if possible */

	/* Parameter adjustments */
	--tau;
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	t_dim1 = *ldt;
	t_offset = 1 + t_dim1;
	t -= t_offset;
	y_dim1 = *ldy;
	y_offset = 1 + y_dim1;
	y -= y_offset;

	/* Function Body */
	if (*n <= 1)
	{
		return 0;
	}

	i__1 = *nb;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		if (i__ > 1)
		{

			/*           Update A(1:n,i) */

			/*           Compute i-th column of A - Y * V' */

			i__2 = i__ - 1;
			dgemv_ ("No transpose", n, &i__2, &c_b347, &y[y_offset], ldy, &a[*k + i__ - 1 + a_dim1], lda, &c_b348, &a[i__ * a_dim1 + 1], &c__1, (ftnlen) 12);

			/*           Apply I - V * T' * V' to this column (call it b) from the */
			/*           left, using the last column of T as workspace */

			/*           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows) */
			/*                    ( V2 )             ( b2 ) */

			/*           where V1 is unit lower triangular */

			/*           w := V1' * b1 */

			i__2 = i__ - 1;
			dcopy_ (&i__2, &a[*k + 1 + i__ * a_dim1], &c__1, &t[*nb * t_dim1 + 1], &c__1);
			i__2 = i__ - 1;
			dtrmv_ ("Lower", "Transpose", "Unit", &i__2, &a[*k + 1 + a_dim1], lda, &t[*nb * t_dim1 + 1], &c__1, (ftnlen) 5, (ftnlen) 9, (ftnlen) 4);

			/*           w := w + V2'*b2 */

			i__2 = *n - *k - i__ + 1;
			i__3 = i__ - 1;
			dgemv_ ("Transpose", &i__2, &i__3, &c_b348, &a[*k + i__ + a_dim1],
					  lda, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b348, &t[*nb * t_dim1 + 1], &c__1, (ftnlen) 9);

			/*           w := T'*w */

			i__2 = i__ - 1;
			dtrmv_ ("Upper", "Transpose", "Non-unit", &i__2, &t[t_offset], ldt, &t[*nb * t_dim1 + 1], &c__1, (ftnlen) 5, (ftnlen) 9, (ftnlen) 8);

			/*           b2 := b2 - V2*w */

			i__2 = *n - *k - i__ + 1;
			i__3 = i__ - 1;
			dgemv_ ("No transpose", &i__2, &i__3, &c_b347, &a[*k + i__ +
																			  a_dim1], lda, &t[*nb * t_dim1 + 1], &c__1, &c_b348, &a[*k + i__ + i__ * a_dim1], &c__1, (ftnlen) 12);

			/*           b1 := b1 - V1*w */

			i__2 = i__ - 1;
			dtrmv_ ("Lower", "No transpose", "Unit", &i__2, &a[*k + 1 + a_dim1], lda, &t[*nb * t_dim1 + 1], &c__1, (ftnlen) 5, (ftnlen) 12, (ftnlen) 4);
			i__2 = i__ - 1;
			daxpy_ (&i__2, &c_b347, &t[*nb * t_dim1 + 1], &c__1, &a[*k + 1 + i__ * a_dim1], &c__1);

			a[*k + i__ - 1 + (i__ - 1) * a_dim1] = ei;
		}

		/*        Generate the elementary reflector H(i) to annihilate */
		/*        A(k+i+1:n,i) */

		i__2 = *n - *k - i__ + 1;
		/* Computing MIN */
		i__3 = *k + i__ + 1;
		dlarfg_ (&i__2, &a[*k + i__ + i__ * a_dim1], &a[min (i__3, *n) + i__ * a_dim1], &c__1, &tau[i__]);
		ei = a[*k + i__ + i__ * a_dim1];
		a[*k + i__ + i__ * a_dim1] = 1.;

		/*        Compute  Y(1:n,i) */

		i__2 = *n - *k - i__ + 1;
		dgemv_ ("No transpose", n, &i__2, &c_b348, &a[(i__ + 1) * a_dim1 + 1],
				  lda, &a[*k + i__ + i__ * a_dim1], &c__1, &c_b507, &y[i__ * y_dim1 + 1], &c__1, (ftnlen) 12);
		i__2 = *n - *k - i__ + 1;
		i__3 = i__ - 1;
		dgemv_ ("Transpose", &i__2, &i__3, &c_b348, &a[*k + i__ + a_dim1], lda,
				  &a[*k + i__ + i__ * a_dim1], &c__1, &c_b507, &t[i__ * t_dim1 + 1], &c__1, (ftnlen) 9);
		i__2 = i__ - 1;
		dgemv_ ("No transpose", n, &i__2, &c_b347, &y[y_offset], ldy, &t[i__ * t_dim1 + 1], &c__1, &c_b348, &y[i__ * y_dim1 + 1], &c__1, (ftnlen) 12);
		dscal_ (n, &tau[i__], &y[i__ * y_dim1 + 1], &c__1);

		/*        Compute T(1:i,i) */

		i__2 = i__ - 1;
		d__1 = -tau[i__];
		dscal_ (&i__2, &d__1, &t[i__ * t_dim1 + 1], &c__1);
		i__2 = i__ - 1;
		dtrmv_ ("Upper", "No transpose", "Non-unit", &i__2, &t[t_offset], ldt, &t[i__ * t_dim1 + 1], &c__1, (ftnlen) 5, (ftnlen) 12, (ftnlen) 8);
		t[i__ + i__ * t_dim1] = tau[i__];

		/* L10: */
	}
	a[*k + *nb + *nb * a_dim1] = ei;

	return 0;

	/*     End of DLAHRD */

}	/* dlahrd_ */

static int
dlaln2_ (logical * ltrans, integer * na, integer * nw,
			doublereal * smin, doublereal * ca, doublereal * a, integer * lda,
			doublereal * d1, doublereal * d2, doublereal * b, integer * ldb,
			doublereal * wr, doublereal * wi, doublereal * x, integer * ldx, doublereal * scale, doublereal * xnorm, integer * info)
{
	/* Initialized data */

	static logical zswap[4] = { FALSE_, FALSE_, TRUE_, TRUE_ };
	static logical rswap[4] = { FALSE_, TRUE_, FALSE_, TRUE_ };
	static integer ipivot[16] /* was [4][4] */  = { 1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2,
		4, 3, 2, 1
	};

	/* System generated locals */
	integer a_dim1, a_offset, b_dim1, b_offset, x_dim1, x_offset;
	doublereal d__1, d__2, d__3, d__4, d__5, d__6;
	static doublereal equiv_0[4], equiv_1[4];

	/* Local variables */
	static integer j;
#define ci (equiv_0)
#define cr (equiv_1)
	static doublereal bi1, bi2, br1, br2, xi1, xi2, xr1, xr2, ci21, ci22, cr21, cr22, li21, csi, ui11, lr21, ui12, ui22;
#define civ (equiv_0)
	static doublereal csr, ur11, ur12, ur22;
#define crv (equiv_1)
	static doublereal bbnd, cmax, ui11r, ui12s, temp, ur11r, ur12s, u22abs;
	static integer icmax;
	static doublereal bnorm, cnorm, smini;
	static doublereal bignum, smlnum;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     October 31, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLALN2 solves a system of the form  (ca A - w D ) X = s B */
	/*  or (ca A' - w D) X = s B   with possible scaling ("s") and */
	/*  perturbation of A.  (A' means A-transpose.) */

	/*  A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA */
	/*  real diagonal matrix, w is a real or complex value, and X and B are */
	/*  NA x 1 matrices -- real if w is real, complex if w is complex.  NA */
	/*  may be 1 or 2. */

	/*  If w is complex, X and B are represented as NA x 2 matrices, */
	/*  the first column of each being the real part and the second */
	/*  being the imaginary part. */

	/*  "s" is a scaling factor (.LE. 1), computed by DLALN2, which is */
	/*  so chosen that X can be computed without overflow.  X is further */
	/*  scaled if necessary to assure that norm(ca A - w D)*norm(X) is less */
	/*  than overflow. */

	/*  If both singular values of (ca A - w D) are less than SMIN, */
	/*  SMIN*identity will be used instead of (ca A - w D).  If only one */
	/*  singular value is less than SMIN, one element of (ca A - w D) will be */
	/*  perturbed enough to make the smallest singular value roughly SMIN. */
	/*  If both singular values are at least SMIN, (ca A - w D) will not be */
	/*  perturbed.  In any case, the perturbation will be at most some small */
	/*  multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values */
	/*  are computed by infinity-norm approximations, and thus will only be */
	/*  correct to a factor of 2 or so. */

	/*  Note: all input quantities are assumed to be smaller than overflow */
	/*  by a reasonable factor.  (See BIGNUM.) */

	/*  Arguments */
	/*  ========== */

	/*  LTRANS  (input) LOGICAL */
	/*          =.TRUE.:  A-transpose will be used. */
	/*          =.FALSE.: A will be used (not transposed.) */

	/*  NA      (input) INTEGER */
	/*          The size of the matrix A.  It may (only) be 1 or 2. */

	/*  NW      (input) INTEGER */
	/*          1 if "w" is real, 2 if "w" is complex.  It may only be 1 */
	/*          or 2. */

	/*  SMIN    (input) DOUBLE PRECISION */
	/*          The desired lower bound on the singular values of A.  This */
	/*          should be a safe distance away from underflow or overflow, */
	/*          say, between (underflow/machine precision) and  (machine */
	/*          precision * overflow ).  (See BIGNUM and ULP.) */

	/*  CA      (input) DOUBLE PRECISION */
	/*          The coefficient c, which A is multiplied by. */

	/*  A       (input) DOUBLE PRECISION array, dimension (LDA,NA) */
	/*          The NA x NA matrix A. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of A.  It must be at least NA. */

	/*  D1      (input) DOUBLE PRECISION */
	/*          The 1,1 element in the diagonal matrix D. */

	/*  D2      (input) DOUBLE PRECISION */
	/*          The 2,2 element in the diagonal matrix D.  Not used if NW=1. */

	/*  B       (input) DOUBLE PRECISION array, dimension (LDB,NW) */
	/*          The NA x NW matrix B (right-hand side).  If NW=2 ("w" is */
	/*          complex), column 1 contains the real part of B and column 2 */
	/*          contains the imaginary part. */

	/*  LDB     (input) INTEGER */
	/*          The leading dimension of B.  It must be at least NA. */

	/*  WR      (input) DOUBLE PRECISION */
	/*          The real part of the scalar "w". */

	/*  WI      (input) DOUBLE PRECISION */
	/*          The imaginary part of the scalar "w".  Not used if NW=1. */

	/*  X       (output) DOUBLE PRECISION array, dimension (LDX,NW) */
	/*          The NA x NW matrix X (unknowns), as computed by DLALN2. */
	/*          If NW=2 ("w" is complex), on exit, column 1 will contain */
	/*          the real part of X and column 2 will contain the imaginary */
	/*          part. */

	/*  LDX     (input) INTEGER */
	/*          The leading dimension of X.  It must be at least NA. */

	/*  SCALE   (output) DOUBLE PRECISION */
	/*          The scale factor that B must be multiplied by to insure */
	/*          that overflow does not occur when computing X.  Thus, */
	/*          (ca A - w D) X  will be SCALE*B, not B (ignoring */
	/*          perturbations of A.)  It will be at most 1. */

	/*  XNORM   (output) DOUBLE PRECISION */
	/*          The infinity-norm of X, when X is regarded as an NA x NW */
	/*          real matrix. */

	/*  INFO    (output) INTEGER */
	/*          An error flag.  It will be set to zero if no error occurs, */
	/*          a negative number if an argument is in error, or a positive */
	/*          number if  ca A - w D  had to be perturbed. */
	/*          The possible values are: */
	/*          = 0: No error occurred, and (ca A - w D) did not have to be */
	/*                 perturbed. */
	/*          = 1: (ca A - w D) had to be perturbed to make its smallest */
	/*               (or only) singular value greater than SMIN. */
	/*          NOTE: In the interests of speed, this routine does not */
	/*                check the inputs for errors. */

	/* ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Local Arrays .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Equivalences .. */
	/*     .. */
	/*     .. Data statements .. */
	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	b_dim1 = *ldb;
	b_offset = 1 + b_dim1;
	b -= b_offset;
	x_dim1 = *ldx;
	x_offset = 1 + x_dim1;
	x -= x_offset;

	/* Function Body */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Compute BIGNUM */

	smlnum = 2. * dlamch_ ("Safe minimum", (ftnlen) 12);
	bignum = 1. / smlnum;
	smini = max (*smin, smlnum);

	/*     Don't check for input errors */

	*info = 0;

	/*     Standard Initializations */

	*scale = 1.;

	if (*na == 1)
	{

		/*        1 x 1  (i.e., scalar) system   C X = B */

		if (*nw == 1)
		{

			/*           Real 1x1 system. */

			/*           C = ca A - w D */

			csr = *ca * a[a_dim1 + 1] - *wr * *d1;
			cnorm = abs (csr);

			/*           If | C | < SMINI, use C = SMINI */

			if (cnorm < smini)
			{
				csr = smini;
				cnorm = smini;
				*info = 1;
			}

			/*           Check scaling for  X = B / C */

			bnorm = (d__1 = b[b_dim1 + 1], abs (d__1));
			if (cnorm < 1. && bnorm > 1.)
			{
				if (bnorm > bignum * cnorm)
				{
					*scale = 1. / bnorm;
				}
			}

			/*           Compute X */

			x[x_dim1 + 1] = b[b_dim1 + 1] * *scale / csr;
			*xnorm = (d__1 = x[x_dim1 + 1], abs (d__1));
		}
		else
		{

			/*           Complex 1x1 system (w is complex) */

			/*           C = ca A - w D */

			csr = *ca * a[a_dim1 + 1] - *wr * *d1;
			csi = -(*wi) * *d1;
			cnorm = abs (csr) + abs (csi);

			/*           If | C | < SMINI, use C = SMINI */

			if (cnorm < smini)
			{
				csr = smini;
				csi = 0.;
				cnorm = smini;
				*info = 1;
			}

			/*           Check scaling for  X = B / C */

			bnorm = (d__1 = b[b_dim1 + 1], abs (d__1)) + (d__2 = b[(b_dim1 << 1) + 1], abs (d__2));
			if (cnorm < 1. && bnorm > 1.)
			{
				if (bnorm > bignum * cnorm)
				{
					*scale = 1. / bnorm;
				}
			}

			/*           Compute X */

			d__1 = *scale * b[b_dim1 + 1];
			d__2 = *scale * b[(b_dim1 << 1) + 1];
			dladiv_ (&d__1, &d__2, &csr, &csi, &x[x_dim1 + 1], &x[(x_dim1 << 1) + 1]);
			*xnorm = (d__1 = x[x_dim1 + 1], abs (d__1)) + (d__2 = x[(x_dim1 << 1) + 1], abs (d__2));
		}

	}
	else
	{

		/*        2x2 System */

		/*        Compute the real part of  C = ca A - w D  (or  ca A' - w D ) */

		cr[0] = *ca * a[a_dim1 + 1] - *wr * *d1;
		cr[3] = *ca * a[(a_dim1 << 1) + 2] - *wr * *d2;
		if (*ltrans)
		{
			cr[2] = *ca * a[a_dim1 + 2];
			cr[1] = *ca * a[(a_dim1 << 1) + 1];
		}
		else
		{
			cr[1] = *ca * a[a_dim1 + 2];
			cr[2] = *ca * a[(a_dim1 << 1) + 1];
		}

		if (*nw == 1)
		{

			/*           Real 2x2 system  (w is real) */

			/*           Find the largest element in C */

			cmax = 0.;
			icmax = 0;

			for (j = 1; j <= 4; ++j)
			{
				if ((d__1 = crv[j - 1], abs (d__1)) > cmax)
				{
					cmax = (d__1 = crv[j - 1], abs (d__1));
					icmax = j;
				}
				/* L10: */
			}

			/*           If norm(C) < SMINI, use SMINI*identity. */

			if (cmax < smini)
			{
				/* Computing MAX */
				d__3 = (d__1 = b[b_dim1 + 1], abs (d__1)), d__4 = (d__2 = b[b_dim1 + 2], abs (d__2));
				bnorm = max (d__3, d__4);
				if (smini < 1. && bnorm > 1.)
				{
					if (bnorm > bignum * smini)
					{
						*scale = 1. / bnorm;
					}
				}
				temp = *scale / smini;
				x[x_dim1 + 1] = temp * b[b_dim1 + 1];
				x[x_dim1 + 2] = temp * b[b_dim1 + 2];
				*xnorm = temp * bnorm;
				*info = 1;
				return 0;
			}

			/*           Gaussian elimination with complete pivoting. */

			ur11 = crv[icmax - 1];
			cr21 = crv[ipivot[(icmax << 2) - 3] - 1];
			ur12 = crv[ipivot[(icmax << 2) - 2] - 1];
			cr22 = crv[ipivot[(icmax << 2) - 1] - 1];
			ur11r = 1. / ur11;
			lr21 = ur11r * cr21;
			ur22 = cr22 - ur12 * lr21;

			/*           If smaller pivot < SMINI, use SMINI */

			if (abs (ur22) < smini)
			{
				ur22 = smini;
				*info = 1;
			}
			if (rswap[icmax - 1])
			{
				br1 = b[b_dim1 + 2];
				br2 = b[b_dim1 + 1];
			}
			else
			{
				br1 = b[b_dim1 + 1];
				br2 = b[b_dim1 + 2];
			}
			br2 -= lr21 * br1;
			/* Computing MAX */
			d__2 = (d__1 = br1 * (ur22 * ur11r), abs (d__1)), d__3 = abs (br2);
			bbnd = max (d__2, d__3);
			if (bbnd > 1. && abs (ur22) < 1.)
			{
				if (bbnd >= bignum * abs (ur22))
				{
					*scale = 1. / bbnd;
				}
			}

			xr2 = br2 * *scale / ur22;
			xr1 = *scale * br1 * ur11r - xr2 * (ur11r * ur12);
			if (zswap[icmax - 1])
			{
				x[x_dim1 + 1] = xr2;
				x[x_dim1 + 2] = xr1;
			}
			else
			{
				x[x_dim1 + 1] = xr1;
				x[x_dim1 + 2] = xr2;
			}
			/* Computing MAX */
			d__1 = abs (xr1), d__2 = abs (xr2);
			*xnorm = max (d__1, d__2);

			/*           Further scaling if  norm(A) norm(X) > overflow */

			if (*xnorm > 1. && cmax > 1.)
			{
				if (*xnorm > bignum / cmax)
				{
					temp = cmax / bignum;
					x[x_dim1 + 1] = temp * x[x_dim1 + 1];
					x[x_dim1 + 2] = temp * x[x_dim1 + 2];
					*xnorm = temp * *xnorm;
					*scale = temp * *scale;
				}
			}
		}
		else
		{

			/*           Complex 2x2 system  (w is complex) */

			/*           Find the largest element in C */

			ci[0] = -(*wi) * *d1;
			ci[1] = 0.;
			ci[2] = 0.;
			ci[3] = -(*wi) * *d2;
			cmax = 0.;
			icmax = 0;

			for (j = 1; j <= 4; ++j)
			{
				if ((d__1 = crv[j - 1], abs (d__1)) + (d__2 = civ[j - 1], abs (d__2)) > cmax)
				{
					cmax = (d__1 = crv[j - 1], abs (d__1)) + (d__2 = civ[j - 1], abs (d__2));
					icmax = j;
				}
				/* L20: */
			}

			/*           If norm(C) < SMINI, use SMINI*identity. */

			if (cmax < smini)
			{
				/* Computing MAX */
				d__5 = (d__1 = b[b_dim1 + 1], abs (d__1)) + (d__2 = b[(b_dim1
																						 << 1) + 1], abs (d__2)), d__6 = (d__3 = b[b_dim1 + 2],
																																	 abs (d__3)) + (d__4 = b[(b_dim1 << 1) + 2], abs (d__4));
				bnorm = max (d__5, d__6);
				if (smini < 1. && bnorm > 1.)
				{
					if (bnorm > bignum * smini)
					{
						*scale = 1. / bnorm;
					}
				}
				temp = *scale / smini;
				x[x_dim1 + 1] = temp * b[b_dim1 + 1];
				x[x_dim1 + 2] = temp * b[b_dim1 + 2];
				x[(x_dim1 << 1) + 1] = temp * b[(b_dim1 << 1) + 1];
				x[(x_dim1 << 1) + 2] = temp * b[(b_dim1 << 1) + 2];
				*xnorm = temp * bnorm;
				*info = 1;
				return 0;
			}

			/*           Gaussian elimination with complete pivoting. */

			ur11 = crv[icmax - 1];
			ui11 = civ[icmax - 1];
			cr21 = crv[ipivot[(icmax << 2) - 3] - 1];
			ci21 = civ[ipivot[(icmax << 2) - 3] - 1];
			ur12 = crv[ipivot[(icmax << 2) - 2] - 1];
			ui12 = civ[ipivot[(icmax << 2) - 2] - 1];
			cr22 = crv[ipivot[(icmax << 2) - 1] - 1];
			ci22 = civ[ipivot[(icmax << 2) - 1] - 1];
			if (icmax == 1 || icmax == 4)
			{

				/*              Code when off-diagonals of pivoted C are real */

				if (abs (ur11) > abs (ui11))
				{
					temp = ui11 / ur11;
					/* Computing 2nd power */
					d__1 = temp;
					ur11r = 1. / (ur11 * (d__1 * d__1 + 1.));
					ui11r = -temp * ur11r;
				}
				else
				{
					temp = ur11 / ui11;
					/* Computing 2nd power */
					d__1 = temp;
					ui11r = -1. / (ui11 * (d__1 * d__1 + 1.));
					ur11r = -temp * ui11r;
				}
				lr21 = cr21 * ur11r;
				li21 = cr21 * ui11r;
				ur12s = ur12 * ur11r;
				ui12s = ur12 * ui11r;
				ur22 = cr22 - ur12 * lr21;
				ui22 = ci22 - ur12 * li21;
			}
			else
			{

				/*              Code when diagonals of pivoted C are real */

				ur11r = 1. / ur11;
				ui11r = 0.;
				lr21 = cr21 * ur11r;
				li21 = ci21 * ur11r;
				ur12s = ur12 * ur11r;
				ui12s = ui12 * ur11r;
				ur22 = cr22 - ur12 * lr21 + ui12 * li21;
				ui22 = -ur12 * li21 - ui12 * lr21;
			}
			u22abs = abs (ur22) + abs (ui22);

			/*           If smaller pivot < SMINI, use SMINI */

			if (u22abs < smini)
			{
				ur22 = smini;
				ui22 = 0.;
				*info = 1;
			}
			if (rswap[icmax - 1])
			{
				br2 = b[b_dim1 + 1];
				br1 = b[b_dim1 + 2];
				bi2 = b[(b_dim1 << 1) + 1];
				bi1 = b[(b_dim1 << 1) + 2];
			}
			else
			{
				br1 = b[b_dim1 + 1];
				br2 = b[b_dim1 + 2];
				bi1 = b[(b_dim1 << 1) + 1];
				bi2 = b[(b_dim1 << 1) + 2];
			}
			br2 = br2 - lr21 * br1 + li21 * bi1;
			bi2 = bi2 - li21 * br1 - lr21 * bi1;
			/* Computing MAX */
			d__1 = (abs (br1) + abs (bi1)) * (u22abs * (abs (ur11r) + abs (ui11r))), d__2 = abs (br2) + abs (bi2);
			bbnd = max (d__1, d__2);
			if (bbnd > 1. && u22abs < 1.)
			{
				if (bbnd >= bignum * u22abs)
				{
					*scale = 1. / bbnd;
					br1 = *scale * br1;
					bi1 = *scale * bi1;
					br2 = *scale * br2;
					bi2 = *scale * bi2;
				}
			}

			dladiv_ (&br2, &bi2, &ur22, &ui22, &xr2, &xi2);
			xr1 = ur11r * br1 - ui11r * bi1 - ur12s * xr2 + ui12s * xi2;
			xi1 = ui11r * br1 + ur11r * bi1 - ui12s * xr2 - ur12s * xi2;
			if (zswap[icmax - 1])
			{
				x[x_dim1 + 1] = xr2;
				x[x_dim1 + 2] = xr1;
				x[(x_dim1 << 1) + 1] = xi2;
				x[(x_dim1 << 1) + 2] = xi1;
			}
			else
			{
				x[x_dim1 + 1] = xr1;
				x[x_dim1 + 2] = xr2;
				x[(x_dim1 << 1) + 1] = xi1;
				x[(x_dim1 << 1) + 2] = xi2;
			}
			/* Computing MAX */
			d__1 = abs (xr1) + abs (xi1), d__2 = abs (xr2) + abs (xi2);
			*xnorm = max (d__1, d__2);

			/*           Further scaling if  norm(A) norm(X) > overflow */

			if (*xnorm > 1. && cmax > 1.)
			{
				if (*xnorm > bignum / cmax)
				{
					temp = cmax / bignum;
					x[x_dim1 + 1] = temp * x[x_dim1 + 1];
					x[x_dim1 + 2] = temp * x[x_dim1 + 2];
					x[(x_dim1 << 1) + 1] = temp * x[(x_dim1 << 1) + 1];
					x[(x_dim1 << 1) + 2] = temp * x[(x_dim1 << 1) + 2];
					*xnorm = temp * *xnorm;
					*scale = temp * *scale;
				}
			}
		}
	}

	return 0;

	/*     End of DLALN2 */

}	/* dlaln2_ */

#undef crv
#undef civ
#undef cr
#undef ci


static doublereal
dlange_ (char *norm, integer * m, integer * n, doublereal * a, integer * lda, doublereal * work, ftnlen norm_len)
{
	/* System generated locals */
	integer a_dim1, a_offset, i__1, i__2;
	doublereal ret_val, d__1, d__2, d__3;

	/* Builtin functions */
	double sqrt (doublereal);

	/* Local variables */
	static integer i__, j;
	static doublereal sum, scale;
	static doublereal value;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     October 31, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLANGE  returns the value of the one norm,  or the Frobenius norm, or */
	/*  the  infinity norm,  or the  element of  largest absolute value  of a */
	/*  real matrix A. */

	/*  Description */
	/*  =========== */

	/*  DLANGE returns the value */

	/*     DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
	/*              ( */
	/*              ( norm1(A),         NORM = '1', 'O' or 'o' */
	/*              ( */
	/*              ( normI(A),         NORM = 'I' or 'i' */
	/*              ( */
	/*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e' */

	/*  where  norm1  denotes the  one norm of a matrix (maximum column sum), */
	/*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and */
	/*  normF  denotes the  Frobenius norm of a matrix (square root of sum of */
	/*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm. */

	/*  Arguments */
	/*  ========= */

	/*  NORM    (input) CHARACTER*1 */
	/*          Specifies the value to be returned in DLANGE as described */
	/*          above. */

	/*  M       (input) INTEGER */
	/*          The number of rows of the matrix A.  M >= 0.  When M = 0, */
	/*          DLANGE is set to zero. */

	/*  N       (input) INTEGER */
	/*          The number of columns of the matrix A.  N >= 0.  When N = 0, */
	/*          DLANGE is set to zero. */

	/*  A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          The m by n matrix A. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(M,1). */

	/*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK), */
	/*          where LWORK >= M when NORM = 'I'; otherwise, WORK is not */
	/*          referenced. */

	/* ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--work;

	/* Function Body */
	if (min (*m, *n) == 0)
	{
		value = 0.;
	}
	else if (lsame_ (norm, "M", (ftnlen) 1, (ftnlen) 1))
	{

		/*        Find max(abs(A(i,j))). */

		value = 0.;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				/* Computing MAX */
				d__2 = value, d__3 = (d__1 = a[i__ + j * a_dim1], abs (d__1));
				value = max (d__2, d__3);
				/* L10: */
			}
			/* L20: */
		}
	}
	else if (lsame_ (norm, "O", (ftnlen) 1, (ftnlen) 1) || *(unsigned char *) norm == '1')
	{

		/*        Find norm1(A). */

		value = 0.;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			sum = 0.;
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				sum += (d__1 = a[i__ + j * a_dim1], abs (d__1));
				/* L30: */
			}
			value = max (value, sum);
			/* L40: */
		}
	}
	else if (lsame_ (norm, "I", (ftnlen) 1, (ftnlen) 1))
	{

		/*        Find normI(A). */

		i__1 = *m;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			work[i__] = 0.;
			/* L50: */
		}
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				work[i__] += (d__1 = a[i__ + j * a_dim1], abs (d__1));
				/* L60: */
			}
			/* L70: */
		}
		value = 0.;
		i__1 = *m;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			/* Computing MAX */
			d__1 = value, d__2 = work[i__];
			value = max (d__1, d__2);
			/* L80: */
		}
	}
	else if (lsame_ (norm, "F", (ftnlen) 1, (ftnlen) 1) || lsame_ (norm, "E", (ftnlen) 1, (ftnlen) 1))
	{

		/*        Find normF(A). */

		scale = 0.;
		sum = 1.;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			dlassq_ (m, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
			/* L90: */
		}
		value = scale * sqrt (sum);
	}

	ret_val = value;
	return ret_val;

	/*     End of DLANGE */

}	/* dlange_ */

static doublereal
dlanhs_ (char *norm, integer * n, doublereal * a, integer * lda, doublereal * work, ftnlen norm_len)
{
	/* System generated locals */
	integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
	doublereal ret_val, d__1, d__2, d__3;

	/* Builtin functions */
	double sqrt (doublereal);

	/* Local variables */
	static integer i__, j;
	static doublereal sum, scale;
	static doublereal value;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     October 31, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLANHS  returns the value of the one norm,  or the Frobenius norm, or */
	/*  the  infinity norm,  or the  element of  largest absolute value  of a */
	/*  Hessenberg matrix A. */

	/*  Description */
	/*  =========== */

	/*  DLANHS returns the value */

	/*     DLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
	/*              ( */
	/*              ( norm1(A),         NORM = '1', 'O' or 'o' */
	/*              ( */
	/*              ( normI(A),         NORM = 'I' or 'i' */
	/*              ( */
	/*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e' */

	/*  where  norm1  denotes the  one norm of a matrix (maximum column sum), */
	/*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and */
	/*  normF  denotes the  Frobenius norm of a matrix (square root of sum of */
	/*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm. */

	/*  Arguments */
	/*  ========= */

	/*  NORM    (input) CHARACTER*1 */
	/*          Specifies the value to be returned in DLANHS as described */
	/*          above. */

	/*  N       (input) INTEGER */
	/*          The order of the matrix A.  N >= 0.  When N = 0, DLANHS is */
	/*          set to zero. */

	/*  A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          The n by n upper Hessenberg matrix A; the part of A below the */
	/*          first sub-diagonal is not referenced. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(N,1). */

	/*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK), */
	/*          where LWORK >= N when NORM = 'I'; otherwise, WORK is not */
	/*          referenced. */

	/* ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--work;

	/* Function Body */
	if (*n == 0)
	{
		value = 0.;
	}
	else if (lsame_ (norm, "M", (ftnlen) 1, (ftnlen) 1))
	{

		/*        Find max(abs(A(i,j))). */

		value = 0.;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			/* Computing MIN */
			i__3 = *n, i__4 = j + 1;
			i__2 = min (i__3, i__4);
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				/* Computing MAX */
				d__2 = value, d__3 = (d__1 = a[i__ + j * a_dim1], abs (d__1));
				value = max (d__2, d__3);
				/* L10: */
			}
			/* L20: */
		}
	}
	else if (lsame_ (norm, "O", (ftnlen) 1, (ftnlen) 1) || *(unsigned char *) norm == '1')
	{

		/*        Find norm1(A). */

		value = 0.;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			sum = 0.;
			/* Computing MIN */
			i__3 = *n, i__4 = j + 1;
			i__2 = min (i__3, i__4);
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				sum += (d__1 = a[i__ + j * a_dim1], abs (d__1));
				/* L30: */
			}
			value = max (value, sum);
			/* L40: */
		}
	}
	else if (lsame_ (norm, "I", (ftnlen) 1, (ftnlen) 1))
	{

		/*        Find normI(A). */

		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			work[i__] = 0.;
			/* L50: */
		}
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			/* Computing MIN */
			i__3 = *n, i__4 = j + 1;
			i__2 = min (i__3, i__4);
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				work[i__] += (d__1 = a[i__ + j * a_dim1], abs (d__1));
				/* L60: */
			}
			/* L70: */
		}
		value = 0.;
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			/* Computing MAX */
			d__1 = value, d__2 = work[i__];
			value = max (d__1, d__2);
			/* L80: */
		}
	}
	else if (lsame_ (norm, "F", (ftnlen) 1, (ftnlen) 1) || lsame_ (norm, "E", (ftnlen) 1, (ftnlen) 1))
	{

		/*        Find normF(A). */

		scale = 0.;
		sum = 1.;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			/* Computing MIN */
			i__3 = *n, i__4 = j + 1;
			i__2 = min (i__3, i__4);
			dlassq_ (&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
			/* L90: */
		}
		value = scale * sqrt (sum);
	}

	ret_val = value;
	return ret_val;

	/*     End of DLANHS */

}	/* dlanhs_ */

static doublereal
dlantr_ (char *norm, char *uplo, char *diag, integer * m, integer * n,
			doublereal * a, integer * lda, doublereal * work, ftnlen norm_len, ftnlen uplo_len, ftnlen diag_len)
{
	/* System generated locals */
	integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
	doublereal ret_val, d__1, d__2, d__3;

	/* Builtin functions */
	double sqrt (doublereal);

	/* Local variables */
	static integer i__, j;
	static doublereal sum, scale;
	static logical udiag;
	static doublereal value;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     October 31, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLANTR  returns the value of the one norm,  or the Frobenius norm, or */
	/*  the  infinity norm,  or the  element of  largest absolute value  of a */
	/*  trapezoidal or triangular matrix A. */

	/*  Description */
	/*  =========== */

	/*  DLANTR returns the value */

	/*     DLANTR = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
	/*              ( */
	/*              ( norm1(A),         NORM = '1', 'O' or 'o' */
	/*              ( */
	/*              ( normI(A),         NORM = 'I' or 'i' */
	/*              ( */
	/*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e' */

	/*  where  norm1  denotes the  one norm of a matrix (maximum column sum), */
	/*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and */
	/*  normF  denotes the  Frobenius norm of a matrix (square root of sum of */
	/*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm. */

	/*  Arguments */
	/*  ========= */

	/*  NORM    (input) CHARACTER*1 */
	/*          Specifies the value to be returned in DLANTR as described */
	/*          above. */

	/*  UPLO    (input) CHARACTER*1 */
	/*          Specifies whether the matrix A is upper or lower trapezoidal. */
	/*          = 'U':  Upper trapezoidal */
	/*          = 'L':  Lower trapezoidal */
	/*          Note that A is triangular instead of trapezoidal if M = N. */

	/*  DIAG    (input) CHARACTER*1 */
	/*          Specifies whether or not the matrix A has unit diagonal. */
	/*          = 'N':  Non-unit diagonal */
	/*          = 'U':  Unit diagonal */

	/*  M       (input) INTEGER */
	/*          The number of rows of the matrix A.  M >= 0, and if */
	/*          UPLO = 'U', M <= N.  When M = 0, DLANTR is set to zero. */

	/*  N       (input) INTEGER */
	/*          The number of columns of the matrix A.  N >= 0, and if */
	/*          UPLO = 'L', N <= M.  When N = 0, DLANTR is set to zero. */

	/*  A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          The trapezoidal matrix A (A is triangular if M = N). */
	/*          If UPLO = 'U', the leading m by n upper trapezoidal part of */
	/*          the array A contains the upper trapezoidal matrix, and the */
	/*          strictly lower triangular part of A is not referenced. */
	/*          If UPLO = 'L', the leading m by n lower trapezoidal part of */
	/*          the array A contains the lower trapezoidal matrix, and the */
	/*          strictly upper triangular part of A is not referenced.  Note */
	/*          that when DIAG = 'U', the diagonal elements of A are not */
	/*          referenced and are assumed to be one. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(M,1). */

	/*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK), */
	/*          where LWORK >= M when NORM = 'I'; otherwise, WORK is not */
	/*          referenced. */

	/* ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--work;

	/* Function Body */
	if (min (*m, *n) == 0)
	{
		value = 0.;
	}
	else if (lsame_ (norm, "M", (ftnlen) 1, (ftnlen) 1))
	{

		/*        Find max(abs(A(i,j))). */

		if (lsame_ (diag, "U", (ftnlen) 1, (ftnlen) 1))
		{
			value = 1.;
			if (lsame_ (uplo, "U", (ftnlen) 1, (ftnlen) 1))
			{
				i__1 = *n;
				for (j = 1; j <= i__1; ++j)
				{
					/* Computing MIN */
					i__3 = *m, i__4 = j - 1;
					i__2 = min (i__3, i__4);
					for (i__ = 1; i__ <= i__2; ++i__)
					{
						/* Computing MAX */
						d__2 = value, d__3 = (d__1 = a[i__ + j * a_dim1], abs (d__1));
						value = max (d__2, d__3);
						/* L10: */
					}
					/* L20: */
				}
			}
			else
			{
				i__1 = *n;
				for (j = 1; j <= i__1; ++j)
				{
					i__2 = *m;
					for (i__ = j + 1; i__ <= i__2; ++i__)
					{
						/* Computing MAX */
						d__2 = value, d__3 = (d__1 = a[i__ + j * a_dim1], abs (d__1));
						value = max (d__2, d__3);
						/* L30: */
					}
					/* L40: */
				}
			}
		}
		else
		{
			value = 0.;
			if (lsame_ (uplo, "U", (ftnlen) 1, (ftnlen) 1))
			{
				i__1 = *n;
				for (j = 1; j <= i__1; ++j)
				{
					i__2 = min (*m, j);
					for (i__ = 1; i__ <= i__2; ++i__)
					{
						/* Computing MAX */
						d__2 = value, d__3 = (d__1 = a[i__ + j * a_dim1], abs (d__1));
						value = max (d__2, d__3);
						/* L50: */
					}
					/* L60: */
				}
			}
			else
			{
				i__1 = *n;
				for (j = 1; j <= i__1; ++j)
				{
					i__2 = *m;
					for (i__ = j; i__ <= i__2; ++i__)
					{
						/* Computing MAX */
						d__2 = value, d__3 = (d__1 = a[i__ + j * a_dim1], abs (d__1));
						value = max (d__2, d__3);
						/* L70: */
					}
					/* L80: */
				}
			}
		}
	}
	else if (lsame_ (norm, "O", (ftnlen) 1, (ftnlen) 1) || *(unsigned char *) norm == '1')
	{

		/*        Find norm1(A). */

		value = 0.;
		udiag = lsame_ (diag, "U", (ftnlen) 1, (ftnlen) 1);
		if (lsame_ (uplo, "U", (ftnlen) 1, (ftnlen) 1))
		{
			i__1 = *n;
			for (j = 1; j <= i__1; ++j)
			{
				if (udiag && j <= *m)
				{
					sum = 1.;
					i__2 = j - 1;
					for (i__ = 1; i__ <= i__2; ++i__)
					{
						sum += (d__1 = a[i__ + j * a_dim1], abs (d__1));
						/* L90: */
					}
				}
				else
				{
					sum = 0.;
					i__2 = min (*m, j);
					for (i__ = 1; i__ <= i__2; ++i__)
					{
						sum += (d__1 = a[i__ + j * a_dim1], abs (d__1));
						/* L100: */
					}
				}
				value = max (value, sum);
				/* L110: */
			}
		}
		else
		{
			i__1 = *n;
			for (j = 1; j <= i__1; ++j)
			{
				if (udiag)
				{
					sum = 1.;
					i__2 = *m;
					for (i__ = j + 1; i__ <= i__2; ++i__)
					{
						sum += (d__1 = a[i__ + j * a_dim1], abs (d__1));
						/* L120: */
					}
				}
				else
				{
					sum = 0.;
					i__2 = *m;
					for (i__ = j; i__ <= i__2; ++i__)
					{
						sum += (d__1 = a[i__ + j * a_dim1], abs (d__1));
						/* L130: */
					}
				}
				value = max (value, sum);
				/* L140: */
			}
		}
	}
	else if (lsame_ (norm, "I", (ftnlen) 1, (ftnlen) 1))
	{

		/*        Find normI(A). */

		if (lsame_ (uplo, "U", (ftnlen) 1, (ftnlen) 1))
		{
			if (lsame_ (diag, "U", (ftnlen) 1, (ftnlen) 1))
			{
				i__1 = *m;
				for (i__ = 1; i__ <= i__1; ++i__)
				{
					work[i__] = 1.;
					/* L150: */
				}
				i__1 = *n;
				for (j = 1; j <= i__1; ++j)
				{
					/* Computing MIN */
					i__3 = *m, i__4 = j - 1;
					i__2 = min (i__3, i__4);
					for (i__ = 1; i__ <= i__2; ++i__)
					{
						work[i__] += (d__1 = a[i__ + j * a_dim1], abs (d__1));
						/* L160: */
					}
					/* L170: */
				}
			}
			else
			{
				i__1 = *m;
				for (i__ = 1; i__ <= i__1; ++i__)
				{
					work[i__] = 0.;
					/* L180: */
				}
				i__1 = *n;
				for (j = 1; j <= i__1; ++j)
				{
					i__2 = min (*m, j);
					for (i__ = 1; i__ <= i__2; ++i__)
					{
						work[i__] += (d__1 = a[i__ + j * a_dim1], abs (d__1));
						/* L190: */
					}
					/* L200: */
				}
			}
		}
		else
		{
			if (lsame_ (diag, "U", (ftnlen) 1, (ftnlen) 1))
			{
				i__1 = *n;
				for (i__ = 1; i__ <= i__1; ++i__)
				{
					work[i__] = 1.;
					/* L210: */
				}
				i__1 = *m;
				for (i__ = *n + 1; i__ <= i__1; ++i__)
				{
					work[i__] = 0.;
					/* L220: */
				}
				i__1 = *n;
				for (j = 1; j <= i__1; ++j)
				{
					i__2 = *m;
					for (i__ = j + 1; i__ <= i__2; ++i__)
					{
						work[i__] += (d__1 = a[i__ + j * a_dim1], abs (d__1));
						/* L230: */
					}
					/* L240: */
				}
			}
			else
			{
				i__1 = *m;
				for (i__ = 1; i__ <= i__1; ++i__)
				{
					work[i__] = 0.;
					/* L250: */
				}
				i__1 = *n;
				for (j = 1; j <= i__1; ++j)
				{
					i__2 = *m;
					for (i__ = j; i__ <= i__2; ++i__)
					{
						work[i__] += (d__1 = a[i__ + j * a_dim1], abs (d__1));
						/* L260: */
					}
					/* L270: */
				}
			}
		}
		value = 0.;
		i__1 = *m;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			/* Computing MAX */
			d__1 = value, d__2 = work[i__];
			value = max (d__1, d__2);
			/* L280: */
		}
	}
	else if (lsame_ (norm, "F", (ftnlen) 1, (ftnlen) 1) || lsame_ (norm, "E", (ftnlen) 1, (ftnlen) 1))
	{

		/*        Find normF(A). */

		if (lsame_ (uplo, "U", (ftnlen) 1, (ftnlen) 1))
		{
			if (lsame_ (diag, "U", (ftnlen) 1, (ftnlen) 1))
			{
				scale = 1.;
				sum = (doublereal) min (*m, *n);
				i__1 = *n;
				for (j = 2; j <= i__1; ++j)
				{
					/* Computing MIN */
					i__3 = *m, i__4 = j - 1;
					i__2 = min (i__3, i__4);
					dlassq_ (&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
					/* L290: */
				}
			}
			else
			{
				scale = 0.;
				sum = 1.;
				i__1 = *n;
				for (j = 1; j <= i__1; ++j)
				{
					i__2 = min (*m, j);
					dlassq_ (&i__2, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
					/* L300: */
				}
			}
		}
		else
		{
			if (lsame_ (diag, "U", (ftnlen) 1, (ftnlen) 1))
			{
				scale = 1.;
				sum = (doublereal) min (*m, *n);
				i__1 = *n;
				for (j = 1; j <= i__1; ++j)
				{
					i__2 = *m - j;
					/* Computing MIN */
					i__3 = *m, i__4 = j + 1;
					dlassq_ (&i__2, &a[min (i__3, i__4) + j * a_dim1], &c__1, &scale, &sum);
					/* L310: */
				}
			}
			else
			{
				scale = 0.;
				sum = 1.;
				i__1 = *n;
				for (j = 1; j <= i__1; ++j)
				{
					i__2 = *m - j + 1;
					dlassq_ (&i__2, &a[j + j * a_dim1], &c__1, &scale, &sum);
					/* L320: */
				}
			}
		}
		value = scale * sqrt (sum);
	}

	ret_val = value;
	return ret_val;

	/*     End of DLANTR */

}	/* dlantr_ */

static int
dlanv2_ (doublereal * a, doublereal * b, doublereal * c__,
			doublereal * d__, doublereal * rt1r, doublereal * rt1i, doublereal * rt2r, doublereal * rt2i, doublereal * cs, doublereal * sn)
{
	/* System generated locals */
	doublereal d__1, d__2;

	/* Builtin functions */
	double d_sign (doublereal *, doublereal *), sqrt (doublereal);

	/* Local variables */
	static doublereal p, z__, aa, bb, cc, dd, cs1, sn1, sab, sac, eps, tau, temp, scale, bcmax, bcmis, sigma;


	/*  -- LAPACK driver routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     June 30, 1999 */

	/*     .. Scalar Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric */
	/*  matrix in standard form: */

	/*       [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ] */
	/*       [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ] */

	/*  where either */
	/*  1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or */
	/*  2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex */
	/*  conjugate eigenvalues. */

	/*  Arguments */
	/*  ========= */

	/*  A       (input/output) DOUBLE PRECISION */
	/*  B       (input/output) DOUBLE PRECISION */
	/*  C       (input/output) DOUBLE PRECISION */
	/*  D       (input/output) DOUBLE PRECISION */
	/*          On entry, the elements of the input matrix. */
	/*          On exit, they are overwritten by the elements of the */
	/*          standardised Schur form. */

	/*  RT1R    (output) DOUBLE PRECISION */
	/*  RT1I    (output) DOUBLE PRECISION */
	/*  RT2R    (output) DOUBLE PRECISION */
	/*  RT2I    (output) DOUBLE PRECISION */
	/*          The real and imaginary parts of the eigenvalues. If the */
	/*          eigenvalues are a complex conjugate pair, RT1I > 0. */

	/*  CS      (output) DOUBLE PRECISION */
	/*  SN      (output) DOUBLE PRECISION */
	/*          Parameters of the rotation matrix. */

	/*  Further Details */
	/*  =============== */

	/*  Modified by V. Sima, Research Institute for Informatics, Bucharest, */
	/*  Romania, to reduce the risk of cancellation errors, */
	/*  when computing real eigenvalues, and to ensure, if possible, that */
	/*  abs(RT1R) >= abs(RT2R). */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	eps = dlamch_ ("P", (ftnlen) 1);
	if (*c__ == 0.)
	{
		*cs = 1.;
		*sn = 0.;
		goto L10;

	}
	else if (*b == 0.)
	{

		/*        Swap rows and columns */

		*cs = 0.;
		*sn = 1.;
		temp = *d__;
		*d__ = *a;
		*a = temp;
		*b = -(*c__);
		*c__ = 0.;
		goto L10;
	}
	else if (*a - *d__ == 0. && d_sign (&c_b348, b) != d_sign (&c_b348, c__))
	{
		*cs = 1.;
		*sn = 0.;
		goto L10;
	}
	else
	{

		temp = *a - *d__;
		p = temp * .5;
		/* Computing MAX */
		d__1 = abs (*b), d__2 = abs (*c__);
		bcmax = max (d__1, d__2);
		/* Computing MIN */
		d__1 = abs (*b), d__2 = abs (*c__);
		bcmis = min (d__1, d__2) * d_sign (&c_b348, b) * d_sign (&c_b348, c__);
		/* Computing MAX */
		d__1 = abs (p);
		scale = max (d__1, bcmax);
		z__ = p / scale * p + bcmax / scale * bcmis;

		/*        If Z is of the order of the machine accuracy, postpone the */
		/*        decision on the nature of eigenvalues */

		if (z__ >= eps * 4.)
		{

			/*           Real eigenvalues. Compute A and D. */

			d__1 = sqrt (scale) * sqrt (z__);
			z__ = p + d_sign (&d__1, &p);
			*a = *d__ + z__;
			*d__ -= bcmax / z__ * bcmis;

			/*           Compute B and the rotation matrix */

			tau = dlapy2_ (c__, &z__);
			*cs = z__ / tau;
			*sn = *c__ / tau;
			*b -= *c__;
			*c__ = 0.;
		}
		else
		{

			/*           Complex eigenvalues, or real (almost) equal eigenvalues. */
			/*           Make diagonal elements equal. */

			sigma = *b + *c__;
			tau = dlapy2_ (&sigma, &temp);
			*cs = sqrt ((abs (sigma) / tau + 1.) * .5);
			*sn = -(p / (tau * *cs)) * d_sign (&c_b348, &sigma);

			/*           Compute [ AA  BB ] = [ A  B ] [ CS -SN ] */
			/*                   [ CC  DD ]   [ C  D ] [ SN  CS ] */

			aa = *a * *cs + *b * *sn;
			bb = -(*a) * *sn + *b * *cs;
			cc = *c__ * *cs + *d__ * *sn;
			dd = -(*c__) * *sn + *d__ * *cs;

			/*           Compute [ A  B ] = [ CS  SN ] [ AA  BB ] */
			/*                   [ C  D ]   [-SN  CS ] [ CC  DD ] */

			*a = aa * *cs + cc * *sn;
			*b = bb * *cs + dd * *sn;
			*c__ = -aa * *sn + cc * *cs;
			*d__ = -bb * *sn + dd * *cs;

			temp = (*a + *d__) * .5;
			*a = temp;
			*d__ = temp;

			if (*c__ != 0.)
			{
				if (*b != 0.)
				{
					if (d_sign (&c_b348, b) == d_sign (&c_b348, c__))
					{

						/*                    Real eigenvalues: reduce to upper triangular form */

						sab = sqrt ((abs (*b)));
						sac = sqrt ((abs (*c__)));
						d__1 = sab * sac;
						p = d_sign (&d__1, c__);
						tau = 1. / sqrt ((d__1 = *b + *c__, abs (d__1)));
						*a = temp + p;
						*d__ = temp - p;
						*b -= *c__;
						*c__ = 0.;
						cs1 = sab * tau;
						sn1 = sac * tau;
						temp = *cs * cs1 - *sn * sn1;
						*sn = *cs * sn1 + *sn * cs1;
						*cs = temp;
					}
				}
				else
				{
					*b = -(*c__);
					*c__ = 0.;
					temp = *cs;
					*cs = -(*sn);
					*sn = temp;
				}
			}
		}

	}

 L10:

	/*     Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I). */

	*rt1r = *a;
	*rt2r = *d__;
	if (*c__ == 0.)
	{
		*rt1i = 0.;
		*rt2i = 0.;
	}
	else
	{
		*rt1i = sqrt ((abs (*b))) * sqrt ((abs (*c__)));
		*rt2i = -(*rt1i);
	}
	return 0;

	/*     End of DLANV2 */

}	/* dlanv2_ */

static doublereal
dlapy2_ (doublereal * x, doublereal * y)
{
	/* System generated locals */
	doublereal ret_val, d__1;

	/* Builtin functions */
	double sqrt (doublereal);

	/* Local variables */
	static doublereal w, z__, xabs, yabs;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     October 31, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary */
	/*  overflow. */

	/*  Arguments */
	/*  ========= */

	/*  X       (input) DOUBLE PRECISION */
	/*  Y       (input) DOUBLE PRECISION */
	/*          X and Y specify the values x and y. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	xabs = abs (*x);
	yabs = abs (*y);
	w = max (xabs, yabs);
	z__ = min (xabs, yabs);
	if (z__ == 0.)
	{
		ret_val = w;
	}
	else
	{
		/* Computing 2nd power */
		d__1 = z__ / w;
		ret_val = w * sqrt (d__1 * d__1 + 1.);
	}
	return ret_val;

	/*     End of DLAPY2 */

}	/* dlapy2_ */

static int
dlaqge_ (integer * m, integer * n, doublereal * a, integer *
			lda, doublereal * r__, doublereal * c__, doublereal * rowcnd, doublereal * colcnd, doublereal * amax, char *equed, ftnlen equed_len)
{
	/* System generated locals */
	integer a_dim1, a_offset, i__1, i__2;

	/* Local variables */
	static integer i__, j;
	static doublereal cj, large, small;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     February 29, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLAQGE equilibrates a general M by N matrix A using the row and */
	/*  scaling factors in the vectors R and C. */

	/*  Arguments */
	/*  ========= */

	/*  M       (input) INTEGER */
	/*          The number of rows of the matrix A.  M >= 0. */

	/*  N       (input) INTEGER */
	/*          The number of columns of the matrix A.  N >= 0. */

	/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          On entry, the M by N matrix A. */
	/*          On exit, the equilibrated matrix.  See EQUED for the form of */
	/*          the equilibrated matrix. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(M,1). */

	/*  R       (input) DOUBLE PRECISION array, dimension (M) */
	/*          The row scale factors for A. */

	/*  C       (input) DOUBLE PRECISION array, dimension (N) */
	/*          The column scale factors for A. */

	/*  ROWCND  (input) DOUBLE PRECISION */
	/*          Ratio of the smallest R(i) to the largest R(i). */

	/*  COLCND  (input) DOUBLE PRECISION */
	/*          Ratio of the smallest C(i) to the largest C(i). */

	/*  AMAX    (input) DOUBLE PRECISION */
	/*          Absolute value of largest matrix entry. */

	/*  EQUED   (output) CHARACTER*1 */
	/*          Specifies the form of equilibration that was done. */
	/*          = 'N':  No equilibration */
	/*          = 'R':  Row equilibration, i.e., A has been premultiplied by */
	/*                  diag(R). */
	/*          = 'C':  Column equilibration, i.e., A has been postmultiplied */
	/*                  by diag(C). */
	/*          = 'B':  Both row and column equilibration, i.e., A has been */
	/*                  replaced by diag(R) * A * diag(C). */

	/*  Internal Parameters */
	/*  =================== */

	/*  THRESH is a threshold value used to decide if row or column scaling */
	/*  should be done based on the ratio of the row or column scaling */
	/*  factors.  If ROWCND < THRESH, row scaling is done, and if */
	/*  COLCND < THRESH, column scaling is done. */

	/*  LARGE and SMALL are threshold values used to decide if row scaling */
	/*  should be done based on the absolute size of the largest matrix */
	/*  element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Quick return if possible */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--r__;
	--c__;

	/* Function Body */
	if (*m <= 0 || *n <= 0)
	{
		*(unsigned char *) equed = 'N';
		return 0;
	}

	/*     Initialize LARGE and SMALL. */

	small = dlamch_ ("Safe minimum", (ftnlen) 12) / dlamch_ ("Precision", (ftnlen) 9);
	large = 1. / small;

	if (*rowcnd >= .1 && *amax >= small && *amax <= large)
	{

		/*        No row scaling */

		if (*colcnd >= .1)
		{

			/*           No column scaling */

			*(unsigned char *) equed = 'N';
		}
		else
		{

			/*           Column scaling */

			i__1 = *n;
			for (j = 1; j <= i__1; ++j)
			{
				cj = c__[j];
				i__2 = *m;
				for (i__ = 1; i__ <= i__2; ++i__)
				{
					a[i__ + j * a_dim1] = cj * a[i__ + j * a_dim1];
					/* L10: */
				}
				/* L20: */
			}
			*(unsigned char *) equed = 'C';
		}
	}
	else if (*colcnd >= .1)
	{

		/*        Row scaling, no column scaling */

		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				a[i__ + j * a_dim1] = r__[i__] * a[i__ + j * a_dim1];
				/* L30: */
			}
			/* L40: */
		}
		*(unsigned char *) equed = 'R';
	}
	else
	{

		/*        Row and column scaling */

		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			cj = c__[j];
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				a[i__ + j * a_dim1] = cj * r__[i__] * a[i__ + j * a_dim1];
				/* L50: */
			}
			/* L60: */
		}
		*(unsigned char *) equed = 'B';
	}

	return 0;

	/*     End of DLAQGE */

}	/* dlaqge_ */

static int
dlaqtr_ (logical * ltran, logical * lreal, integer * n,
			doublereal * t, integer * ldt, doublereal * b, doublereal * w, doublereal * scale, doublereal * x, doublereal * work, integer * info)
{
	/* System generated locals */
	integer t_dim1, t_offset, i__1, i__2;
	doublereal d__1, d__2, d__3, d__4, d__5, d__6;

	/* Local variables */
	static doublereal d__[4] /* was [2][2] */ ;
	static integer i__, j, k;
	static doublereal v[4] /* was [2][2] */ , z__;
	static integer j1, j2, n1, n2;
	static doublereal si, xj, sr, rec, eps, tjj, tmp;
	static integer ierr;
	static doublereal smin, xmax;
	static integer jnext;
	static doublereal sminw, xnorm;
	static doublereal scaloc;
	static doublereal bignum;
	static logical notran;
	static doublereal smlnum;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     June 30, 1999 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLAQTR solves the real quasi-triangular system */

	/*               op(T)*p = scale*c,               if LREAL = .TRUE. */

	/*  or the complex quasi-triangular systems */

	/*             op(T + iB)*(p+iq) = scale*(c+id),  if LREAL = .FALSE. */

	/*  in real arithmetic, where T is upper quasi-triangular. */
	/*  If LREAL = .FALSE., then the first diagonal block of T must be */
	/*  1 by 1, B is the specially structured matrix */

	/*                 B = [ b(1) b(2) ... b(n) ] */
	/*                     [       w            ] */
	/*                     [           w        ] */
	/*                     [              .     ] */
	/*                     [                 w  ] */

	/*  op(A) = A or A', A' denotes the conjugate transpose of */
	/*  matrix A. */

	/*  On input, X = [ c ].  On output, X = [ p ]. */
	/*                [ d ]                  [ q ] */

	/*  This subroutine is designed for the condition number estimation */
	/*  in routine DTRSNA. */

	/*  Arguments */
	/*  ========= */

	/*  LTRAN   (input) LOGICAL */
	/*          On entry, LTRAN specifies the option of conjugate transpose: */
	/*             = .FALSE.,    op(T+i*B) = T+i*B, */
	/*             = .TRUE.,     op(T+i*B) = (T+i*B)'. */

	/*  LREAL   (input) LOGICAL */
	/*          On entry, LREAL specifies the input matrix structure: */
	/*             = .FALSE.,    the input is complex */
	/*             = .TRUE.,     the input is real */

	/*  N       (input) INTEGER */
	/*          On entry, N specifies the order of T+i*B. N >= 0. */

	/*  T       (input) DOUBLE PRECISION array, dimension (LDT,N) */
	/*          On entry, T contains a matrix in Schur canonical form. */
	/*          If LREAL = .FALSE., then the first diagonal block of T mu */
	/*          be 1 by 1. */

	/*  LDT     (input) INTEGER */
	/*          The leading dimension of the matrix T. LDT >= max(1,N). */

	/*  B       (input) DOUBLE PRECISION array, dimension (N) */
	/*          On entry, B contains the elements to form the matrix */
	/*          B as described above. */
	/*          If LREAL = .TRUE., B is not referenced. */

	/*  W       (input) DOUBLE PRECISION */
	/*          On entry, W is the diagonal element of the matrix B. */
	/*          If LREAL = .TRUE., W is not referenced. */

	/*  SCALE   (output) DOUBLE PRECISION */
	/*          On exit, SCALE is the scale factor. */

	/*  X       (input/output) DOUBLE PRECISION array, dimension (2*N) */
	/*          On entry, X contains the right hand side of the system. */
	/*          On exit, X is overwritten by the solution. */

	/*  WORK    (workspace) DOUBLE PRECISION array, dimension (N) */

	/*  INFO    (output) INTEGER */
	/*          On exit, INFO is set to */
	/*             0: successful exit. */
	/*               1: the some diagonal 1 by 1 block has been perturbed by */
	/*                  a small number SMIN to keep nonsingularity. */
	/*               2: the some diagonal 2 by 2 block has been perturbed by */
	/*                  a small number in DLALN2 to keep nonsingularity. */
	/*          NOTE: In the interests of speed, this routine does not */
	/*                check the inputs for errors. */

	/* ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Local Arrays .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Do not test the input parameters for errors */

	/* Parameter adjustments */
	t_dim1 = *ldt;
	t_offset = 1 + t_dim1;
	t -= t_offset;
	--b;
	--x;
	--work;

	/* Function Body */
	notran = !(*ltran);
	*info = 0;

	/*     Quick return if possible */

	if (*n == 0)
	{
		return 0;
	}

	/*     Set constants to control overflow */

	eps = dlamch_ ("P", (ftnlen) 1);
	smlnum = dlamch_ ("S", (ftnlen) 1) / eps;
	bignum = 1. / smlnum;

	xnorm = dlange_ ("M", n, n, &t[t_offset], ldt, d__, (ftnlen) 1);
	if (!(*lreal))
	{
		/* Computing MAX */
		d__1 = xnorm, d__2 = abs (*w), d__1 = max (d__1, d__2), d__2 = dlange_ ("M", n, &c__1, &b[1], n, d__, (ftnlen) 1);
		xnorm = max (d__1, d__2);
	}
	/* Computing MAX */
	d__1 = smlnum, d__2 = eps * xnorm;
	smin = max (d__1, d__2);

	/*     Compute 1-norm of each column of strictly upper triangular */
	/*     part of T to control overflow in triangular solver. */

	work[1] = 0.;
	i__1 = *n;
	for (j = 2; j <= i__1; ++j)
	{
		i__2 = j - 1;
		work[j] = dasum_ (&i__2, &t[j * t_dim1 + 1], &c__1);
		/* L10: */
	}

	if (!(*lreal))
	{
		i__1 = *n;
		for (i__ = 2; i__ <= i__1; ++i__)
		{
			work[i__] += (d__1 = b[i__], abs (d__1));
			/* L20: */
		}
	}

	n2 = *n << 1;
	n1 = *n;
	if (!(*lreal))
	{
		n1 = n2;
	}
	k = idamax_ (&n1, &x[1], &c__1);
	xmax = (d__1 = x[k], abs (d__1));
	*scale = 1.;

	if (xmax > bignum)
	{
		*scale = bignum / xmax;
		dscal_ (&n1, scale, &x[1], &c__1);
		xmax = bignum;
	}

	if (*lreal)
	{

		if (notran)
		{

			/*           Solve T*p = scale*c */

			jnext = *n;
			for (j = *n; j >= 1; --j)
			{
				if (j > jnext)
				{
					goto L30;
				}
				j1 = j;
				j2 = j;
				jnext = j - 1;
				if (j > 1)
				{
					if (t[j + (j - 1) * t_dim1] != 0.)
					{
						j1 = j - 1;
						jnext = j - 2;
					}
				}

				if (j1 == j2)
				{

					/*                 Meet 1 by 1 diagonal block */

					/*                 Scale to avoid overflow when computing */
					/*                     x(j) = b(j)/T(j,j) */

					xj = (d__1 = x[j1], abs (d__1));
					tjj = (d__1 = t[j1 + j1 * t_dim1], abs (d__1));
					tmp = t[j1 + j1 * t_dim1];
					if (tjj < smin)
					{
						tmp = smin;
						tjj = smin;
						*info = 1;
					}

					if (xj == 0.)
					{
						goto L30;
					}

					if (tjj < 1.)
					{
						if (xj > bignum * tjj)
						{
							rec = 1. / xj;
							dscal_ (n, &rec, &x[1], &c__1);
							*scale *= rec;
							xmax *= rec;
						}
					}
					x[j1] /= tmp;
					xj = (d__1 = x[j1], abs (d__1));

					/*                 Scale x if necessary to avoid overflow when adding a */
					/*                 multiple of column j1 of T. */

					if (xj > 1.)
					{
						rec = 1. / xj;
						if (work[j1] > (bignum - xmax) * rec)
						{
							dscal_ (n, &rec, &x[1], &c__1);
							*scale *= rec;
						}
					}
					if (j1 > 1)
					{
						i__1 = j1 - 1;
						d__1 = -x[j1];
						daxpy_ (&i__1, &d__1, &t[j1 * t_dim1 + 1], &c__1, &x[1], &c__1);
						i__1 = j1 - 1;
						k = idamax_ (&i__1, &x[1], &c__1);
						xmax = (d__1 = x[k], abs (d__1));
					}

				}
				else
				{

					/*                 Meet 2 by 2 diagonal block */

					/*                 Call 2 by 2 linear system solve, to take */
					/*                 care of possible overflow by scaling factor. */

					d__[0] = x[j1];
					d__[1] = x[j2];
					dlaln2_ (&c_false, &c__2, &c__1, &smin, &c_b348, &t[j1 +
																						 j1 * t_dim1], ldt, &c_b348, &c_b348, d__, &c__2, &c_b507, &c_b507, v, &c__2, &scaloc, &xnorm,
								&ierr);
					if (ierr != 0)
					{
						*info = 2;
					}

					if (scaloc != 1.)
					{
						dscal_ (n, &scaloc, &x[1], &c__1);
						*scale *= scaloc;
					}
					x[j1] = v[0];
					x[j2] = v[1];

					/*                 Scale V(1,1) (= X(J1)) and/or V(2,1) (=X(J2)) */
					/*                 to avoid overflow in updating right-hand side. */

					/* Computing MAX */
					d__1 = abs (v[0]), d__2 = abs (v[1]);
					xj = max (d__1, d__2);
					if (xj > 1.)
					{
						rec = 1. / xj;
						/* Computing MAX */
						d__1 = work[j1], d__2 = work[j2];
						if (max (d__1, d__2) > (bignum - xmax) * rec)
						{
							dscal_ (n, &rec, &x[1], &c__1);
							*scale *= rec;
						}
					}

					/*                 Update right-hand side */

					if (j1 > 1)
					{
						i__1 = j1 - 1;
						d__1 = -x[j1];
						daxpy_ (&i__1, &d__1, &t[j1 * t_dim1 + 1], &c__1, &x[1], &c__1);
						i__1 = j1 - 1;
						d__1 = -x[j2];
						daxpy_ (&i__1, &d__1, &t[j2 * t_dim1 + 1], &c__1, &x[1], &c__1);
						i__1 = j1 - 1;
						k = idamax_ (&i__1, &x[1], &c__1);
						xmax = (d__1 = x[k], abs (d__1));
					}

				}

			 L30:
				;
			}

		}
		else
		{

			/*           Solve T'*p = scale*c */

			jnext = 1;
			i__1 = *n;
			for (j = 1; j <= i__1; ++j)
			{
				if (j < jnext)
				{
					goto L40;
				}
				j1 = j;
				j2 = j;
				jnext = j + 1;
				if (j < *n)
				{
					if (t[j + 1 + j * t_dim1] != 0.)
					{
						j2 = j + 1;
						jnext = j + 2;
					}
				}

				if (j1 == j2)
				{

					/*                 1 by 1 diagonal block */

					/*                 Scale if necessary to avoid overflow in forming the */
					/*                 right-hand side element by inner product. */

					xj = (d__1 = x[j1], abs (d__1));
					if (xmax > 1.)
					{
						rec = 1. / xmax;
						if (work[j1] > (bignum - xj) * rec)
						{
							dscal_ (n, &rec, &x[1], &c__1);
							*scale *= rec;
							xmax *= rec;
						}
					}

					i__2 = j1 - 1;
					x[j1] -= ddot_ (&i__2, &t[j1 * t_dim1 + 1], &c__1, &x[1], &c__1);

					xj = (d__1 = x[j1], abs (d__1));
					tjj = (d__1 = t[j1 + j1 * t_dim1], abs (d__1));
					tmp = t[j1 + j1 * t_dim1];
					if (tjj < smin)
					{
						tmp = smin;
						tjj = smin;
						*info = 1;
					}

					if (tjj < 1.)
					{
						if (xj > bignum * tjj)
						{
							rec = 1. / xj;
							dscal_ (n, &rec, &x[1], &c__1);
							*scale *= rec;
							xmax *= rec;
						}
					}
					x[j1] /= tmp;
					/* Computing MAX */
					d__2 = xmax, d__3 = (d__1 = x[j1], abs (d__1));
					xmax = max (d__2, d__3);

				}
				else
				{

					/*                 2 by 2 diagonal block */

					/*                 Scale if necessary to avoid overflow in forming the */
					/*                 right-hand side elements by inner product. */

					/* Computing MAX */
					d__3 = (d__1 = x[j1], abs (d__1)), d__4 = (d__2 = x[j2], abs (d__2));
					xj = max (d__3, d__4);
					if (xmax > 1.)
					{
						rec = 1. / xmax;
						/* Computing MAX */
						d__1 = work[j2], d__2 = work[j1];
						if (max (d__1, d__2) > (bignum - xj) * rec)
						{
							dscal_ (n, &rec, &x[1], &c__1);
							*scale *= rec;
							xmax *= rec;
						}
					}

					i__2 = j1 - 1;
					d__[0] = x[j1] - ddot_ (&i__2, &t[j1 * t_dim1 + 1], &c__1, &x[1], &c__1);
					i__2 = j1 - 1;
					d__[1] = x[j2] - ddot_ (&i__2, &t[j2 * t_dim1 + 1], &c__1, &x[1], &c__1);

					dlaln2_ (&c_true, &c__2, &c__1, &smin, &c_b348, &t[j1 + j1
																						* t_dim1], ldt, &c_b348, &c_b348, d__, &c__2, &c_b507, &c_b507, v, &c__2, &scaloc, &xnorm,
								&ierr);
					if (ierr != 0)
					{
						*info = 2;
					}

					if (scaloc != 1.)
					{
						dscal_ (n, &scaloc, &x[1], &c__1);
						*scale *= scaloc;
					}
					x[j1] = v[0];
					x[j2] = v[1];
					/* Computing MAX */
					d__3 = (d__1 = x[j1], abs (d__1)), d__4 = (d__2 = x[j2], abs (d__2)), d__3 = max (d__3, d__4);
					xmax = max (d__3, xmax);

				}
			 L40:
				;
			}
		}

	}
	else
	{

		/* Computing MAX */
		d__1 = eps * abs (*w);
		sminw = max (d__1, smin);
		if (notran)
		{

			/*           Solve (T + iB)*(p+iq) = c+id */

			jnext = *n;
			for (j = *n; j >= 1; --j)
			{
				if (j > jnext)
				{
					goto L70;
				}
				j1 = j;
				j2 = j;
				jnext = j - 1;
				if (j > 1)
				{
					if (t[j + (j - 1) * t_dim1] != 0.)
					{
						j1 = j - 1;
						jnext = j - 2;
					}
				}

				if (j1 == j2)
				{

					/*                 1 by 1 diagonal block */

					/*                 Scale if necessary to avoid overflow in division */

					z__ = *w;
					if (j1 == 1)
					{
						z__ = b[1];
					}
					xj = (d__1 = x[j1], abs (d__1)) + (d__2 = x[*n + j1], abs (d__2));
					tjj = (d__1 = t[j1 + j1 * t_dim1], abs (d__1)) + abs (z__);
					tmp = t[j1 + j1 * t_dim1];
					if (tjj < sminw)
					{
						tmp = sminw;
						tjj = sminw;
						*info = 1;
					}

					if (xj == 0.)
					{
						goto L70;
					}

					if (tjj < 1.)
					{
						if (xj > bignum * tjj)
						{
							rec = 1. / xj;
							dscal_ (&n2, &rec, &x[1], &c__1);
							*scale *= rec;
							xmax *= rec;
						}
					}
					dladiv_ (&x[j1], &x[*n + j1], &tmp, &z__, &sr, &si);
					x[j1] = sr;
					x[*n + j1] = si;
					xj = (d__1 = x[j1], abs (d__1)) + (d__2 = x[*n + j1], abs (d__2));

					/*                 Scale x if necessary to avoid overflow when adding a */
					/*                 multiple of column j1 of T. */

					if (xj > 1.)
					{
						rec = 1. / xj;
						if (work[j1] > (bignum - xmax) * rec)
						{
							dscal_ (&n2, &rec, &x[1], &c__1);
							*scale *= rec;
						}
					}

					if (j1 > 1)
					{
						i__1 = j1 - 1;
						d__1 = -x[j1];
						daxpy_ (&i__1, &d__1, &t[j1 * t_dim1 + 1], &c__1, &x[1], &c__1);
						i__1 = j1 - 1;
						d__1 = -x[*n + j1];
						daxpy_ (&i__1, &d__1, &t[j1 * t_dim1 + 1], &c__1, &x[*n + 1], &c__1);

						x[1] += b[j1] * x[*n + j1];
						x[*n + 1] -= b[j1] * x[j1];

						xmax = 0.;
						i__1 = j1 - 1;
						for (k = 1; k <= i__1; ++k)
						{
							/* Computing MAX */
							d__3 = xmax, d__4 = (d__1 = x[k], abs (d__1)) + (d__2 = x[k + *n], abs (d__2));
							xmax = max (d__3, d__4);
							/* L50: */
						}
					}

				}
				else
				{

					/*                 Meet 2 by 2 diagonal block */

					d__[0] = x[j1];
					d__[1] = x[j2];
					d__[2] = x[*n + j1];
					d__[3] = x[*n + j2];
					d__1 = -(*w);
					dlaln2_ (&c_false, &c__2, &c__2, &sminw, &c_b348, &t[j1 +
																						  j1 * t_dim1], ldt, &c_b348, &c_b348, d__, &c__2, &c_b507, &d__1, v, &c__2, &scaloc, &xnorm,
								&ierr);
					if (ierr != 0)
					{
						*info = 2;
					}

					if (scaloc != 1.)
					{
						i__1 = *n << 1;
						dscal_ (&i__1, &scaloc, &x[1], &c__1);
						*scale = scaloc * *scale;
					}
					x[j1] = v[0];
					x[j2] = v[1];
					x[*n + j1] = v[2];
					x[*n + j2] = v[3];

					/*                 Scale X(J1), .... to avoid overflow in */
					/*                 updating right hand side. */

					/* Computing MAX */
					d__1 = abs (v[0]) + abs (v[2]), d__2 = abs (v[1]) + abs (v[3]);
					xj = max (d__1, d__2);
					if (xj > 1.)
					{
						rec = 1. / xj;
						/* Computing MAX */
						d__1 = work[j1], d__2 = work[j2];
						if (max (d__1, d__2) > (bignum - xmax) * rec)
						{
							dscal_ (&n2, &rec, &x[1], &c__1);
							*scale *= rec;
						}
					}

					/*                 Update the right-hand side. */

					if (j1 > 1)
					{
						i__1 = j1 - 1;
						d__1 = -x[j1];
						daxpy_ (&i__1, &d__1, &t[j1 * t_dim1 + 1], &c__1, &x[1], &c__1);
						i__1 = j1 - 1;
						d__1 = -x[j2];
						daxpy_ (&i__1, &d__1, &t[j2 * t_dim1 + 1], &c__1, &x[1], &c__1);

						i__1 = j1 - 1;
						d__1 = -x[*n + j1];
						daxpy_ (&i__1, &d__1, &t[j1 * t_dim1 + 1], &c__1, &x[*n + 1], &c__1);
						i__1 = j1 - 1;
						d__1 = -x[*n + j2];
						daxpy_ (&i__1, &d__1, &t[j2 * t_dim1 + 1], &c__1, &x[*n + 1], &c__1);

						x[1] = x[1] + b[j1] * x[*n + j1] + b[j2] * x[*n + j2];
						x[*n + 1] = x[*n + 1] - b[j1] * x[j1] - b[j2] * x[j2];

						xmax = 0.;
						i__1 = j1 - 1;
						for (k = 1; k <= i__1; ++k)
						{
							/* Computing MAX */
							d__3 = (d__1 = x[k], abs (d__1)) + (d__2 = x[k + *n], abs (d__2));
							xmax = max (d__3, xmax);
							/* L60: */
						}
					}

				}
			 L70:
				;
			}

		}
		else
		{

			/*           Solve (T + iB)'*(p+iq) = c+id */

			jnext = 1;
			i__1 = *n;
			for (j = 1; j <= i__1; ++j)
			{
				if (j < jnext)
				{
					goto L80;
				}
				j1 = j;
				j2 = j;
				jnext = j + 1;
				if (j < *n)
				{
					if (t[j + 1 + j * t_dim1] != 0.)
					{
						j2 = j + 1;
						jnext = j + 2;
					}
				}

				if (j1 == j2)
				{

					/*                 1 by 1 diagonal block */

					/*                 Scale if necessary to avoid overflow in forming the */
					/*                 right-hand side element by inner product. */

					xj = (d__1 = x[j1], abs (d__1)) + (d__2 = x[j1 + *n], abs (d__2));
					if (xmax > 1.)
					{
						rec = 1. / xmax;
						if (work[j1] > (bignum - xj) * rec)
						{
							dscal_ (&n2, &rec, &x[1], &c__1);
							*scale *= rec;
							xmax *= rec;
						}
					}

					i__2 = j1 - 1;
					x[j1] -= ddot_ (&i__2, &t[j1 * t_dim1 + 1], &c__1, &x[1], &c__1);
					i__2 = j1 - 1;
					x[*n + j1] -= ddot_ (&i__2, &t[j1 * t_dim1 + 1], &c__1, &x[*n + 1], &c__1);
					if (j1 > 1)
					{
						x[j1] -= b[j1] * x[*n + 1];
						x[*n + j1] += b[j1] * x[1];
					}
					xj = (d__1 = x[j1], abs (d__1)) + (d__2 = x[j1 + *n], abs (d__2));

					z__ = *w;
					if (j1 == 1)
					{
						z__ = b[1];
					}

					/*                 Scale if necessary to avoid overflow in */
					/*                 complex division */

					tjj = (d__1 = t[j1 + j1 * t_dim1], abs (d__1)) + abs (z__);
					tmp = t[j1 + j1 * t_dim1];
					if (tjj < sminw)
					{
						tmp = sminw;
						tjj = sminw;
						*info = 1;
					}

					if (tjj < 1.)
					{
						if (xj > bignum * tjj)
						{
							rec = 1. / xj;
							dscal_ (&n2, &rec, &x[1], &c__1);
							*scale *= rec;
							xmax *= rec;
						}
					}
					d__1 = -z__;
					dladiv_ (&x[j1], &x[*n + j1], &tmp, &d__1, &sr, &si);
					x[j1] = sr;
					x[j1 + *n] = si;
					/* Computing MAX */
					d__3 = (d__1 = x[j1], abs (d__1)) + (d__2 = x[j1 + *n], abs (d__2));
					xmax = max (d__3, xmax);

				}
				else
				{

					/*                 2 by 2 diagonal block */

					/*                 Scale if necessary to avoid overflow in forming the */
					/*                 right-hand side element by inner product. */

					/* Computing MAX */
					d__5 = (d__1 = x[j1], abs (d__1)) + (d__2 = x[*n + j1], abs (d__2)), d__6 = (d__3 = x[j2], abs (d__3)) + (d__4 = x[*n + j2], abs (d__4));
					xj = max (d__5, d__6);
					if (xmax > 1.)
					{
						rec = 1. / xmax;
						/* Computing MAX */
						d__1 = work[j1], d__2 = work[j2];
						if (max (d__1, d__2) > (bignum - xj) / xmax)
						{
							dscal_ (&n2, &rec, &x[1], &c__1);
							*scale *= rec;
							xmax *= rec;
						}
					}

					i__2 = j1 - 1;
					d__[0] = x[j1] - ddot_ (&i__2, &t[j1 * t_dim1 + 1], &c__1, &x[1], &c__1);
					i__2 = j1 - 1;
					d__[1] = x[j2] - ddot_ (&i__2, &t[j2 * t_dim1 + 1], &c__1, &x[1], &c__1);
					i__2 = j1 - 1;
					d__[2] = x[*n + j1] - ddot_ (&i__2, &t[j1 * t_dim1 + 1], &c__1, &x[*n + 1], &c__1);
					i__2 = j1 - 1;
					d__[3] = x[*n + j2] - ddot_ (&i__2, &t[j2 * t_dim1 + 1], &c__1, &x[*n + 1], &c__1);
					d__[0] -= b[j1] * x[*n + 1];
					d__[1] -= b[j2] * x[*n + 1];
					d__[2] += b[j1] * x[1];
					d__[3] += b[j2] * x[1];

					dlaln2_ (&c_true, &c__2, &c__2, &sminw, &c_b348, &t[j1 +
																						 j1 * t_dim1], ldt, &c_b348, &c_b348, d__, &c__2, &c_b507, w, v, &c__2, &scaloc, &xnorm,
								&ierr);
					if (ierr != 0)
					{
						*info = 2;
					}

					if (scaloc != 1.)
					{
						dscal_ (&n2, &scaloc, &x[1], &c__1);
						*scale = scaloc * *scale;
					}
					x[j1] = v[0];
					x[j2] = v[1];
					x[*n + j1] = v[2];
					x[*n + j2] = v[3];
					/* Computing MAX */
					d__5 = (d__1 = x[j1], abs (d__1)) + (d__2 = x[*n + j1],
																	 abs (d__2)), d__6 = (d__3 = x[j2], abs (d__3)) + (d__4 = x[*n + j2], abs (d__4)), d__5 = max (d__5, d__6);
					xmax = max (d__5, xmax);

				}

			 L80:
				;
			}

		}

	}

	return 0;

	/*     End of DLAQTR */

}	/* dlaqtr_ */

static int
dlarfb_ (char *side, char *trans, char *direct, char *storev, integer * m, integer * n, integer * k, doublereal * v, integer *
			ldv, doublereal * t, integer * ldt, doublereal * c__, integer * ldc,
			doublereal * work, integer * ldwork, ftnlen side_len, ftnlen trans_len, ftnlen direct_len, ftnlen storev_len)
{
	/* System generated locals */
	integer c_dim1, c_offset, t_dim1, t_offset, v_dim1, v_offset, work_dim1, work_offset, i__1, i__2;

	/* Local variables */
	static integer i__, j;
	static char transt[1];


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     February 29, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLARFB applies a real block reflector H or its transpose H' to a */
	/*  real m by n matrix C, from either the left or the right. */

	/*  Arguments */
	/*  ========= */

	/*  SIDE    (input) CHARACTER*1 */
	/*          = 'L': apply H or H' from the Left */
	/*          = 'R': apply H or H' from the Right */

	/*  TRANS   (input) CHARACTER*1 */
	/*          = 'N': apply H (No transpose) */
	/*          = 'T': apply H' (Transpose) */

	/*  DIRECT  (input) CHARACTER*1 */
	/*          Indicates how H is formed from a product of elementary */
	/*          reflectors */
	/*          = 'F': H = H(1) H(2) . . . H(k) (Forward) */
	/*          = 'B': H = H(k) . . . H(2) H(1) (Backward) */

	/*  STOREV  (input) CHARACTER*1 */
	/*          Indicates how the vectors which define the elementary */
	/*          reflectors are stored: */
	/*          = 'C': Columnwise */
	/*          = 'R': Rowwise */

	/*  M       (input) INTEGER */
	/*          The number of rows of the matrix C. */

	/*  N       (input) INTEGER */
	/*          The number of columns of the matrix C. */

	/*  K       (input) INTEGER */
	/*          The order of the matrix T (= the number of elementary */
	/*          reflectors whose product defines the block reflector). */

	/*  V       (input) DOUBLE PRECISION array, dimension */
	/*                                (LDV,K) if STOREV = 'C' */
	/*                                (LDV,M) if STOREV = 'R' and SIDE = 'L' */
	/*                                (LDV,N) if STOREV = 'R' and SIDE = 'R' */
	/*          The matrix V. See further details. */

	/*  LDV     (input) INTEGER */
	/*          The leading dimension of the array V. */
	/*          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M); */
	/*          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N); */
	/*          if STOREV = 'R', LDV >= K. */

	/*  T       (input) DOUBLE PRECISION array, dimension (LDT,K) */
	/*          The triangular k by k matrix T in the representation of the */
	/*          block reflector. */

	/*  LDT     (input) INTEGER */
	/*          The leading dimension of the array T. LDT >= K. */

	/*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
	/*          On entry, the m by n matrix C. */
	/*          On exit, C is overwritten by H*C or H'*C or C*H or C*H'. */

	/*  LDC     (input) INTEGER */
	/*          The leading dimension of the array C. LDA >= max(1,M). */

	/*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,K) */

	/*  LDWORK  (input) INTEGER */
	/*          The leading dimension of the array WORK. */
	/*          If SIDE = 'L', LDWORK >= max(1,N); */
	/*          if SIDE = 'R', LDWORK >= max(1,M). */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Quick return if possible */

	/* Parameter adjustments */
	v_dim1 = *ldv;
	v_offset = 1 + v_dim1;
	v -= v_offset;
	t_dim1 = *ldt;
	t_offset = 1 + t_dim1;
	t -= t_offset;
	c_dim1 = *ldc;
	c_offset = 1 + c_dim1;
	c__ -= c_offset;
	work_dim1 = *ldwork;
	work_offset = 1 + work_dim1;
	work -= work_offset;

	/* Function Body */
	if (*m <= 0 || *n <= 0)
	{
		return 0;
	}

	if (lsame_ (trans, "N", (ftnlen) 1, (ftnlen) 1))
	{
		*(unsigned char *) transt = 'T';
	}
	else
	{
		*(unsigned char *) transt = 'N';
	}

	if (lsame_ (storev, "C", (ftnlen) 1, (ftnlen) 1))
	{

		if (lsame_ (direct, "F", (ftnlen) 1, (ftnlen) 1))
		{

			/*           Let  V =  ( V1 )    (first K rows) */
			/*                     ( V2 ) */
			/*           where  V1  is unit lower triangular. */

			if (lsame_ (side, "L", (ftnlen) 1, (ftnlen) 1))
			{

				/*              Form  H * C  or  H' * C  where  C = ( C1 ) */
				/*                                                  ( C2 ) */

				/*              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK) */

				/*              W := C1' */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j)
				{
					dcopy_ (n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1], &c__1);
					/* L10: */
				}

				/*              W := W * V1 */

				dtrmm_ ("Right", "Lower", "No transpose", "Unit", n, k, &c_b348, &v[v_offset], ldv, &work[work_offset], ldwork,
						  (ftnlen) 5, (ftnlen) 5, (ftnlen) 12, (ftnlen) 4);
				if (*m > *k)
				{

					/*                 W := W + C2'*V2 */

					i__1 = *m - *k;
					dgemm_ ("Transpose", "No transpose", n, k, &i__1, &c_b348,
							  &c__[*k + 1 + c_dim1], ldc, &v[*k + 1 + v_dim1], ldv, &c_b348, &work[work_offset], ldwork, (ftnlen) 9, (ftnlen) 12);
				}

				/*              W := W * T'  or  W * T */

				dtrmm_ ("Right", "Upper", transt, "Non-unit", n, k, &c_b348, &t[t_offset], ldt, &work[work_offset], ldwork, (ftnlen)
						  5, (ftnlen) 5, (ftnlen) 1, (ftnlen) 8);

				/*              C := C - V * W' */

				if (*m > *k)
				{

					/*                 C2 := C2 - V2 * W' */

					i__1 = *m - *k;
					dgemm_ ("No transpose", "Transpose", &i__1, n, k, &c_b347,
							  &v[*k + 1 + v_dim1], ldv, &work[work_offset], ldwork, &c_b348, &c__[*k + 1 + c_dim1], ldc, (ftnlen) 12, (ftnlen) 9);
				}

				/*              W := W * V1' */

				dtrmm_ ("Right", "Lower", "Transpose", "Unit", n, k, &c_b348, &v[v_offset], ldv, &work[work_offset], ldwork, (ftnlen)
						  5, (ftnlen) 5, (ftnlen) 9, (ftnlen) 4);

				/*              C1 := C1 - W' */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j)
				{
					i__2 = *n;
					for (i__ = 1; i__ <= i__2; ++i__)
					{
						c__[j + i__ * c_dim1] -= work[i__ + j * work_dim1];
						/* L20: */
					}
					/* L30: */
				}

			}
			else if (lsame_ (side, "R", (ftnlen) 1, (ftnlen) 1))
			{

				/*              Form  C * H  or  C * H'  where  C = ( C1  C2 ) */

				/*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK) */

				/*              W := C1 */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j)
				{
					dcopy_ (m, &c__[j * c_dim1 + 1], &c__1, &work[j * work_dim1 + 1], &c__1);
					/* L40: */
				}

				/*              W := W * V1 */

				dtrmm_ ("Right", "Lower", "No transpose", "Unit", m, k, &c_b348, &v[v_offset], ldv, &work[work_offset], ldwork,
						  (ftnlen) 5, (ftnlen) 5, (ftnlen) 12, (ftnlen) 4);
				if (*n > *k)
				{

					/*                 W := W + C2 * V2 */

					i__1 = *n - *k;
					dgemm_ ("No transpose", "No transpose", m, k, &i__1, &c_b348, &c__[(*k + 1) * c_dim1 + 1], ldc, &v[*k +
																																						1 + v_dim1], ldv, &c_b348, &work[work_offset],
							  ldwork, (ftnlen) 12, (ftnlen) 12);
				}

				/*              W := W * T  or  W * T' */

				dtrmm_ ("Right", "Upper", trans, "Non-unit", m, k, &c_b348, &t[t_offset], ldt, &work[work_offset], ldwork, (ftnlen) 5,
						  (ftnlen) 5, (ftnlen) 1, (ftnlen) 8);

				/*              C := C - W * V' */

				if (*n > *k)
				{

					/*                 C2 := C2 - W * V2' */

					i__1 = *n - *k;
					dgemm_ ("No transpose", "Transpose", m, &i__1, k, &c_b347,
							  &work[work_offset], ldwork, &v[*k + 1 + v_dim1], ldv, &c_b348, &c__[(*k + 1) * c_dim1 + 1], ldc, (ftnlen) 12, (ftnlen) 9);
				}

				/*              W := W * V1' */

				dtrmm_ ("Right", "Lower", "Transpose", "Unit", m, k, &c_b348, &v[v_offset], ldv, &work[work_offset], ldwork, (ftnlen)
						  5, (ftnlen) 5, (ftnlen) 9, (ftnlen) 4);

				/*              C1 := C1 - W */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j)
				{
					i__2 = *m;
					for (i__ = 1; i__ <= i__2; ++i__)
					{
						c__[i__ + j * c_dim1] -= work[i__ + j * work_dim1];
						/* L50: */
					}
					/* L60: */
				}
			}

		}
		else
		{

			/*           Let  V =  ( V1 ) */
			/*                     ( V2 )    (last K rows) */
			/*           where  V2  is unit upper triangular. */

			if (lsame_ (side, "L", (ftnlen) 1, (ftnlen) 1))
			{

				/*              Form  H * C  or  H' * C  where  C = ( C1 ) */
				/*                                                  ( C2 ) */

				/*              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK) */

				/*              W := C2' */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j)
				{
					dcopy_ (n, &c__[*m - *k + j + c_dim1], ldc, &work[j * work_dim1 + 1], &c__1);
					/* L70: */
				}

				/*              W := W * V2 */

				dtrmm_ ("Right", "Upper", "No transpose", "Unit", n, k, &c_b348, &v[*m - *k + 1 + v_dim1], ldv, &work[work_offset], ldwork, (ftnlen) 5, (ftnlen) 5,
						  (ftnlen) 12, (ftnlen) 4);
				if (*m > *k)
				{

					/*                 W := W + C1'*V1 */

					i__1 = *m - *k;
					dgemm_ ("Transpose", "No transpose", n, k, &i__1, &c_b348,
							  &c__[c_offset], ldc, &v[v_offset], ldv, &c_b348, &work[work_offset], ldwork, (ftnlen) 9, (ftnlen) 12);
				}

				/*              W := W * T'  or  W * T */

				dtrmm_ ("Right", "Lower", transt, "Non-unit", n, k, &c_b348, &t[t_offset], ldt, &work[work_offset], ldwork, (ftnlen)
						  5, (ftnlen) 5, (ftnlen) 1, (ftnlen) 8);

				/*              C := C - V * W' */

				if (*m > *k)
				{

					/*                 C1 := C1 - V1 * W' */

					i__1 = *m - *k;
					dgemm_ ("No transpose", "Transpose", &i__1, n, k, &c_b347,
							  &v[v_offset], ldv, &work[work_offset], ldwork, &c_b348, &c__[c_offset], ldc, (ftnlen) 12, (ftnlen) 9);
				}

				/*              W := W * V2' */

				dtrmm_ ("Right", "Upper", "Transpose", "Unit", n, k, &c_b348, &v[*m - *k + 1 + v_dim1], ldv, &work[work_offset],
						  ldwork, (ftnlen) 5, (ftnlen) 5, (ftnlen) 9, (ftnlen) 4);

				/*              C2 := C2 - W' */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j)
				{
					i__2 = *n;
					for (i__ = 1; i__ <= i__2; ++i__)
					{
						c__[*m - *k + j + i__ * c_dim1] -= work[i__ + j * work_dim1];
						/* L80: */
					}
					/* L90: */
				}

			}
			else if (lsame_ (side, "R", (ftnlen) 1, (ftnlen) 1))
			{

				/*              Form  C * H  or  C * H'  where  C = ( C1  C2 ) */

				/*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK) */

				/*              W := C2 */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j)
				{
					dcopy_ (m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[j * work_dim1 + 1], &c__1);
					/* L100: */
				}

				/*              W := W * V2 */

				dtrmm_ ("Right", "Upper", "No transpose", "Unit", m, k, &c_b348, &v[*n - *k + 1 + v_dim1], ldv, &work[work_offset], ldwork, (ftnlen) 5, (ftnlen) 5,
						  (ftnlen) 12, (ftnlen) 4);
				if (*n > *k)
				{

					/*                 W := W + C1 * V1 */

					i__1 = *n - *k;
					dgemm_ ("No transpose", "No transpose", m, k, &i__1, &c_b348, &c__[c_offset], ldc, &v[v_offset], ldv, &c_b348, &work[work_offset], ldwork,
							  (ftnlen) 12, (ftnlen) 12);
				}

				/*              W := W * T  or  W * T' */

				dtrmm_ ("Right", "Lower", trans, "Non-unit", m, k, &c_b348, &t[t_offset], ldt, &work[work_offset], ldwork, (ftnlen) 5,
						  (ftnlen) 5, (ftnlen) 1, (ftnlen) 8);

				/*              C := C - W * V' */

				if (*n > *k)
				{

					/*                 C1 := C1 - W * V1' */

					i__1 = *n - *k;
					dgemm_ ("No transpose", "Transpose", m, &i__1, k, &c_b347,
							  &work[work_offset], ldwork, &v[v_offset], ldv, &c_b348, &c__[c_offset], ldc, (ftnlen) 12, (ftnlen) 9);
				}

				/*              W := W * V2' */

				dtrmm_ ("Right", "Upper", "Transpose", "Unit", m, k, &c_b348, &v[*n - *k + 1 + v_dim1], ldv, &work[work_offset],
						  ldwork, (ftnlen) 5, (ftnlen) 5, (ftnlen) 9, (ftnlen) 4);

				/*              C2 := C2 - W */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j)
				{
					i__2 = *m;
					for (i__ = 1; i__ <= i__2; ++i__)
					{
						c__[i__ + (*n - *k + j) * c_dim1] -= work[i__ + j * work_dim1];
						/* L110: */
					}
					/* L120: */
				}
			}
		}

	}
	else if (lsame_ (storev, "R", (ftnlen) 1, (ftnlen) 1))
	{

		if (lsame_ (direct, "F", (ftnlen) 1, (ftnlen) 1))
		{

			/*           Let  V =  ( V1  V2 )    (V1: first K columns) */
			/*           where  V1  is unit upper triangular. */

			if (lsame_ (side, "L", (ftnlen) 1, (ftnlen) 1))
			{

				/*              Form  H * C  or  H' * C  where  C = ( C1 ) */
				/*                                                  ( C2 ) */

				/*              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK) */

				/*              W := C1' */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j)
				{
					dcopy_ (n, &c__[j + c_dim1], ldc, &work[j * work_dim1 + 1], &c__1);
					/* L130: */
				}

				/*              W := W * V1' */

				dtrmm_ ("Right", "Upper", "Transpose", "Unit", n, k, &c_b348, &v[v_offset], ldv, &work[work_offset], ldwork, (ftnlen)
						  5, (ftnlen) 5, (ftnlen) 9, (ftnlen) 4);
				if (*m > *k)
				{

					/*                 W := W + C2'*V2' */

					i__1 = *m - *k;
					dgemm_ ("Transpose", "Transpose", n, k, &i__1, &c_b348, &c__[*k + 1 + c_dim1], ldc, &v[(*k + 1) * v_dim1 +
																																		1], ldv, &c_b348, &work[work_offset], ldwork, (ftnlen) 9,
							  (ftnlen) 9);
				}

				/*              W := W * T'  or  W * T */

				dtrmm_ ("Right", "Upper", transt, "Non-unit", n, k, &c_b348, &t[t_offset], ldt, &work[work_offset], ldwork, (ftnlen)
						  5, (ftnlen) 5, (ftnlen) 1, (ftnlen) 8);

				/*              C := C - V' * W' */

				if (*m > *k)
				{

					/*                 C2 := C2 - V2' * W' */

					i__1 = *m - *k;
					dgemm_ ("Transpose", "Transpose", &i__1, n, k, &c_b347, &v[(*k + 1) * v_dim1 + 1], ldv, &work[work_offset],
							  ldwork, &c_b348, &c__[*k + 1 + c_dim1], ldc, (ftnlen) 9, (ftnlen) 9);
				}

				/*              W := W * V1 */

				dtrmm_ ("Right", "Upper", "No transpose", "Unit", n, k, &c_b348, &v[v_offset], ldv, &work[work_offset], ldwork,
						  (ftnlen) 5, (ftnlen) 5, (ftnlen) 12, (ftnlen) 4);

				/*              C1 := C1 - W' */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j)
				{
					i__2 = *n;
					for (i__ = 1; i__ <= i__2; ++i__)
					{
						c__[j + i__ * c_dim1] -= work[i__ + j * work_dim1];
						/* L140: */
					}
					/* L150: */
				}

			}
			else if (lsame_ (side, "R", (ftnlen) 1, (ftnlen) 1))
			{

				/*              Form  C * H  or  C * H'  where  C = ( C1  C2 ) */

				/*              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK) */

				/*              W := C1 */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j)
				{
					dcopy_ (m, &c__[j * c_dim1 + 1], &c__1, &work[j * work_dim1 + 1], &c__1);
					/* L160: */
				}

				/*              W := W * V1' */

				dtrmm_ ("Right", "Upper", "Transpose", "Unit", m, k, &c_b348, &v[v_offset], ldv, &work[work_offset], ldwork, (ftnlen)
						  5, (ftnlen) 5, (ftnlen) 9, (ftnlen) 4);
				if (*n > *k)
				{

					/*                 W := W + C2 * V2' */

					i__1 = *n - *k;
					dgemm_ ("No transpose", "Transpose", m, k, &i__1, &c_b348,
							  &c__[(*k + 1) * c_dim1 + 1], ldc, &v[(*k + 1) * v_dim1 + 1], ldv, &c_b348, &work[work_offset], ldwork, (ftnlen) 12, (ftnlen) 9);
				}

				/*              W := W * T  or  W * T' */

				dtrmm_ ("Right", "Upper", trans, "Non-unit", m, k, &c_b348, &t[t_offset], ldt, &work[work_offset], ldwork, (ftnlen) 5,
						  (ftnlen) 5, (ftnlen) 1, (ftnlen) 8);

				/*              C := C - W * V */

				if (*n > *k)
				{

					/*                 C2 := C2 - W * V2 */

					i__1 = *n - *k;
					dgemm_ ("No transpose", "No transpose", m, &i__1, k, &c_b347, &work[work_offset], ldwork, &v[(*k + 1) *
																																				v_dim1 + 1], ldv, &c_b348, &c__[(*k + 1) * c_dim1
																																														  + 1], ldc,
							  (ftnlen) 12, (ftnlen) 12);
				}

				/*              W := W * V1 */

				dtrmm_ ("Right", "Upper", "No transpose", "Unit", m, k, &c_b348, &v[v_offset], ldv, &work[work_offset], ldwork,
						  (ftnlen) 5, (ftnlen) 5, (ftnlen) 12, (ftnlen) 4);

				/*              C1 := C1 - W */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j)
				{
					i__2 = *m;
					for (i__ = 1; i__ <= i__2; ++i__)
					{
						c__[i__ + j * c_dim1] -= work[i__ + j * work_dim1];
						/* L170: */
					}
					/* L180: */
				}

			}

		}
		else
		{

			/*           Let  V =  ( V1  V2 )    (V2: last K columns) */
			/*           where  V2  is unit lower triangular. */

			if (lsame_ (side, "L", (ftnlen) 1, (ftnlen) 1))
			{

				/*              Form  H * C  or  H' * C  where  C = ( C1 ) */
				/*                                                  ( C2 ) */

				/*              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK) */

				/*              W := C2' */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j)
				{
					dcopy_ (n, &c__[*m - *k + j + c_dim1], ldc, &work[j * work_dim1 + 1], &c__1);
					/* L190: */
				}

				/*              W := W * V2' */

				dtrmm_ ("Right", "Lower", "Transpose", "Unit", n, k, &c_b348, &v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[work_offset], ldwork, (ftnlen) 5,
						  (ftnlen) 5, (ftnlen) 9, (ftnlen) 4);
				if (*m > *k)
				{

					/*                 W := W + C1'*V1' */

					i__1 = *m - *k;
					dgemm_ ("Transpose", "Transpose", n, k, &i__1, &c_b348, &c__[c_offset], ldc, &v[v_offset], ldv, &c_b348, &work[work_offset], ldwork, (ftnlen) 9,
							  (ftnlen) 9);
				}

				/*              W := W * T'  or  W * T */

				dtrmm_ ("Right", "Lower", transt, "Non-unit", n, k, &c_b348, &t[t_offset], ldt, &work[work_offset], ldwork, (ftnlen)
						  5, (ftnlen) 5, (ftnlen) 1, (ftnlen) 8);

				/*              C := C - V' * W' */

				if (*m > *k)
				{

					/*                 C1 := C1 - V1' * W' */

					i__1 = *m - *k;
					dgemm_ ("Transpose", "Transpose", &i__1, n, k, &c_b347, &v[v_offset], ldv, &work[work_offset], ldwork, &c_b348, &c__[c_offset], ldc, (ftnlen) 9,
							  (ftnlen) 9);
				}

				/*              W := W * V2 */

				dtrmm_ ("Right", "Lower", "No transpose", "Unit", n, k, &c_b348, &v[(*m - *k + 1) * v_dim1 + 1], ldv, &work[work_offset], ldwork, (ftnlen) 5,
						  (ftnlen) 5, (ftnlen) 12, (ftnlen) 4);

				/*              C2 := C2 - W' */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j)
				{
					i__2 = *n;
					for (i__ = 1; i__ <= i__2; ++i__)
					{
						c__[*m - *k + j + i__ * c_dim1] -= work[i__ + j * work_dim1];
						/* L200: */
					}
					/* L210: */
				}

			}
			else if (lsame_ (side, "R", (ftnlen) 1, (ftnlen) 1))
			{

				/*              Form  C * H  or  C * H'  where  C = ( C1  C2 ) */

				/*              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK) */

				/*              W := C2 */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j)
				{
					dcopy_ (m, &c__[(*n - *k + j) * c_dim1 + 1], &c__1, &work[j * work_dim1 + 1], &c__1);
					/* L220: */
				}

				/*              W := W * V2' */

				dtrmm_ ("Right", "Lower", "Transpose", "Unit", m, k, &c_b348, &v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[work_offset], ldwork, (ftnlen) 5,
						  (ftnlen) 5, (ftnlen) 9, (ftnlen) 4);
				if (*n > *k)
				{

					/*                 W := W + C1 * V1' */

					i__1 = *n - *k;
					dgemm_ ("No transpose", "Transpose", m, k, &i__1, &c_b348,
							  &c__[c_offset], ldc, &v[v_offset], ldv, &c_b348, &work[work_offset], ldwork, (ftnlen) 12, (ftnlen) 9);
				}

				/*              W := W * T  or  W * T' */

				dtrmm_ ("Right", "Lower", trans, "Non-unit", m, k, &c_b348, &t[t_offset], ldt, &work[work_offset], ldwork, (ftnlen) 5,
						  (ftnlen) 5, (ftnlen) 1, (ftnlen) 8);

				/*              C := C - W * V */

				if (*n > *k)
				{

					/*                 C1 := C1 - W * V1 */

					i__1 = *n - *k;
					dgemm_ ("No transpose", "No transpose", m, &i__1, k, &c_b347, &work[work_offset], ldwork, &v[v_offset],
							  ldv, &c_b348, &c__[c_offset], ldc, (ftnlen) 12, (ftnlen) 12);
				}

				/*              W := W * V2 */

				dtrmm_ ("Right", "Lower", "No transpose", "Unit", m, k, &c_b348, &v[(*n - *k + 1) * v_dim1 + 1], ldv, &work[work_offset], ldwork, (ftnlen) 5,
						  (ftnlen) 5, (ftnlen) 12, (ftnlen) 4);

				/*              C1 := C1 - W */

				i__1 = *k;
				for (j = 1; j <= i__1; ++j)
				{
					i__2 = *m;
					for (i__ = 1; i__ <= i__2; ++i__)
					{
						c__[i__ + (*n - *k + j) * c_dim1] -= work[i__ + j * work_dim1];
						/* L230: */
					}
					/* L240: */
				}

			}

		}
	}

	return 0;

	/*     End of DLARFB */

}	/* dlarfb_ */

static int
dlarf_ (char *side, integer * m, integer * n, doublereal * v,
		  integer * incv, doublereal * tau, doublereal * c__, integer * ldc, doublereal * work, ftnlen side_len)
{
	/* System generated locals */
	integer c_dim1, c_offset;
	doublereal d__1;

	/* Local variables */


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     February 29, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLARF applies a real elementary reflector H to a real m by n matrix */
	/*  C, from either the left or the right. H is represented in the form */

	/*        H = I - tau * v * v' */

	/*  where tau is a real scalar and v is a real vector. */

	/*  If tau = 0, then H is taken to be the unit matrix. */

	/*  Arguments */
	/*  ========= */

	/*  SIDE    (input) CHARACTER*1 */
	/*          = 'L': form  H * C */
	/*          = 'R': form  C * H */

	/*  M       (input) INTEGER */
	/*          The number of rows of the matrix C. */

	/*  N       (input) INTEGER */
	/*          The number of columns of the matrix C. */

	/*  V       (input) DOUBLE PRECISION array, dimension */
	/*                     (1 + (M-1)*abs(INCV)) if SIDE = 'L' */
	/*                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R' */
	/*          The vector v in the representation of H. V is not used if */
	/*          TAU = 0. */

	/*  INCV    (input) INTEGER */
	/*          The increment between elements of v. INCV <> 0. */

	/*  TAU     (input) DOUBLE PRECISION */
	/*          The value tau in the representation of H. */

	/*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
	/*          On entry, the m by n matrix C. */
	/*          On exit, C is overwritten by the matrix H * C if SIDE = 'L', */
	/*          or C * H if SIDE = 'R'. */

	/*  LDC     (input) INTEGER */
	/*          The leading dimension of the array C. LDC >= max(1,M). */

	/*  WORK    (workspace) DOUBLE PRECISION array, dimension */
	/*                         (N) if SIDE = 'L' */
	/*                      or (M) if SIDE = 'R' */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/* Parameter adjustments */
	--v;
	c_dim1 = *ldc;
	c_offset = 1 + c_dim1;
	c__ -= c_offset;
	--work;

	/* Function Body */
	if (lsame_ (side, "L", (ftnlen) 1, (ftnlen) 1))
	{

		/*        Form  H * C */

		if (*tau != 0.)
		{

			/*           w := C' * v */

			dgemv_ ("Transpose", m, n, &c_b348, &c__[c_offset], ldc, &v[1], incv, &c_b507, &work[1], &c__1, (ftnlen) 9);

			/*           C := C - v * w' */

			d__1 = -(*tau);
			dger_ (m, n, &d__1, &v[1], incv, &work[1], &c__1, &c__[c_offset], ldc);
		}
	}
	else
	{

		/*        Form  C * H */

		if (*tau != 0.)
		{

			/*           w := C * v */

			dgemv_ ("No transpose", m, n, &c_b348, &c__[c_offset], ldc, &v[1], incv, &c_b507, &work[1], &c__1, (ftnlen) 12);

			/*           C := C - w * v' */

			d__1 = -(*tau);
			dger_ (m, n, &d__1, &work[1], &c__1, &v[1], incv, &c__[c_offset], ldc);
		}
	}
	return 0;

	/*     End of DLARF */

}	/* dlarf_ */

static int
dlarfg_ (integer * n, doublereal * alpha, doublereal * x, integer * incx, doublereal * tau)
{
	/* System generated locals */
	integer i__1;
	doublereal d__1;

	/* Builtin functions */
	double d_sign (doublereal *, doublereal *);

	/* Local variables */
	static integer j, knt;
	static doublereal beta;
	static doublereal xnorm;
	static doublereal safmin, rsafmn;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     September 30, 1994 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLARFG generates a real elementary reflector H of order n, such */
	/*  that */

	/*        H * ( alpha ) = ( beta ),   H' * H = I. */
	/*            (   x   )   (   0  ) */

	/*  where alpha and beta are scalars, and x is an (n-1)-element real */
	/*  vector. H is represented in the form */

	/*        H = I - tau * ( 1 ) * ( 1 v' ) , */
	/*                      ( v ) */

	/*  where tau is a real scalar and v is a real (n-1)-element */
	/*  vector. */

	/*  If the elements of x are all zero, then tau = 0 and H is taken to be */
	/*  the unit matrix. */

	/*  Otherwise  1 <= tau <= 2. */

	/*  Arguments */
	/*  ========= */

	/*  N       (input) INTEGER */
	/*          The order of the elementary reflector. */

	/*  ALPHA   (input/output) DOUBLE PRECISION */
	/*          On entry, the value alpha. */
	/*          On exit, it is overwritten with the value beta. */

	/*  X       (input/output) DOUBLE PRECISION array, dimension */
	/*                         (1+(N-2)*abs(INCX)) */
	/*          On entry, the vector x. */
	/*          On exit, it is overwritten with the vector v. */

	/*  INCX    (input) INTEGER */
	/*          The increment between elements of X. INCX > 0. */

	/*  TAU     (output) DOUBLE PRECISION */
	/*          The value tau. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/* Parameter adjustments */
	--x;

	/* Function Body */
	if (*n <= 1)
	{
		*tau = 0.;
		return 0;
	}

	i__1 = *n - 1;
	xnorm = dnrm2_ (&i__1, &x[1], incx);

	if (xnorm == 0.)
	{

		/*        H  =  I */

		*tau = 0.;
	}
	else
	{

		/*        general case */

		d__1 = dlapy2_ (alpha, &xnorm);
		beta = -d_sign (&d__1, alpha);
		safmin = dlamch_ ("S", (ftnlen) 1) / dlamch_ ("E", (ftnlen) 1);
		if (abs (beta) < safmin)
		{

			/*           XNORM, BETA may be inaccurate; scale X and recompute them */

			rsafmn = 1. / safmin;
			knt = 0;
		 L10:
			++knt;
			i__1 = *n - 1;
			dscal_ (&i__1, &rsafmn, &x[1], incx);
			beta *= rsafmn;
			*alpha *= rsafmn;
			if (abs (beta) < safmin)
			{
				goto L10;
			}

			/*           New BETA is at most 1, at least SAFMIN */

			i__1 = *n - 1;
			xnorm = dnrm2_ (&i__1, &x[1], incx);
			d__1 = dlapy2_ (alpha, &xnorm);
			beta = -d_sign (&d__1, alpha);
			*tau = (beta - *alpha) / beta;
			i__1 = *n - 1;
			d__1 = 1. / (*alpha - beta);
			dscal_ (&i__1, &d__1, &x[1], incx);

			/*           If ALPHA is subnormal, it may lose relative accuracy */

			*alpha = beta;
			i__1 = knt;
			for (j = 1; j <= i__1; ++j)
			{
				*alpha *= safmin;
				/* L20: */
			}
		}
		else
		{
			*tau = (beta - *alpha) / beta;
			i__1 = *n - 1;
			d__1 = 1. / (*alpha - beta);
			dscal_ (&i__1, &d__1, &x[1], incx);
			*alpha = beta;
		}
	}

	return 0;

	/*     End of DLARFG */

}	/* dlarfg_ */

static int
dlarft_ (char *direct, char *storev, integer * n, integer *
			k, doublereal * v, integer * ldv, doublereal * tau, doublereal * t, integer * ldt, ftnlen direct_len, ftnlen storev_len)
{
	/* System generated locals */
	integer t_dim1, t_offset, v_dim1, v_offset, i__1, i__2, i__3;
	doublereal d__1;

	/* Local variables */
	static integer i__, j;
	static doublereal vii;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     February 29, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLARFT forms the triangular factor T of a real block reflector H */
	/*  of order n, which is defined as a product of k elementary reflectors. */

	/*  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular; */

	/*  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular. */

	/*  If STOREV = 'C', the vector which defines the elementary reflector */
	/*  H(i) is stored in the i-th column of the array V, and */

	/*     H  =  I - V * T * V' */

	/*  If STOREV = 'R', the vector which defines the elementary reflector */
	/*  H(i) is stored in the i-th row of the array V, and */

	/*     H  =  I - V' * T * V */

	/*  Arguments */
	/*  ========= */

	/*  DIRECT  (input) CHARACTER*1 */
	/*          Specifies the order in which the elementary reflectors are */
	/*          multiplied to form the block reflector: */
	/*          = 'F': H = H(1) H(2) . . . H(k) (Forward) */
	/*          = 'B': H = H(k) . . . H(2) H(1) (Backward) */

	/*  STOREV  (input) CHARACTER*1 */
	/*          Specifies how the vectors which define the elementary */
	/*          reflectors are stored (see also Further Details): */
	/*          = 'C': columnwise */
	/*          = 'R': rowwise */

	/*  N       (input) INTEGER */
	/*          The order of the block reflector H. N >= 0. */

	/*  K       (input) INTEGER */
	/*          The order of the triangular factor T (= the number of */
	/*          elementary reflectors). K >= 1. */

	/*  V       (input/output) DOUBLE PRECISION array, dimension */
	/*                               (LDV,K) if STOREV = 'C' */
	/*                               (LDV,N) if STOREV = 'R' */
	/*          The matrix V. See further details. */

	/*  LDV     (input) INTEGER */
	/*          The leading dimension of the array V. */
	/*          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K. */

	/*  TAU     (input) DOUBLE PRECISION array, dimension (K) */
	/*          TAU(i) must contain the scalar factor of the elementary */
	/*          reflector H(i). */

	/*  T       (output) DOUBLE PRECISION array, dimension (LDT,K) */
	/*          The k by k triangular factor T of the block reflector. */
	/*          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is */
	/*          lower triangular. The rest of the array is not used. */

	/*  LDT     (input) INTEGER */
	/*          The leading dimension of the array T. LDT >= K. */

	/*  Further Details */
	/*  =============== */

	/*  The shape of the matrix V and the storage of the vectors which define */
	/*  the H(i) is best illustrated by the following example with n = 5 and */
	/*  k = 3. The elements equal to 1 are not stored; the corresponding */
	/*  array elements are modified but restored on exit. The rest of the */
	/*  array is not used. */

	/*  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R': */

	/*               V = (  1       )                 V = (  1 v1 v1 v1 v1 ) */
	/*                   ( v1  1    )                     (     1 v2 v2 v2 ) */
	/*                   ( v1 v2  1 )                     (        1 v3 v3 ) */
	/*                   ( v1 v2 v3 ) */
	/*                   ( v1 v2 v3 ) */

	/*  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R': */

	/*               V = ( v1 v2 v3 )                 V = ( v1 v1  1       ) */
	/*                   ( v1 v2 v3 )                     ( v2 v2 v2  1    ) */
	/*                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 ) */
	/*                   (     1 v3 ) */
	/*                   (        1 ) */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Quick return if possible */

	/* Parameter adjustments */
	v_dim1 = *ldv;
	v_offset = 1 + v_dim1;
	v -= v_offset;
	--tau;
	t_dim1 = *ldt;
	t_offset = 1 + t_dim1;
	t -= t_offset;

	/* Function Body */
	if (*n == 0)
	{
		return 0;
	}

	if (lsame_ (direct, "F", (ftnlen) 1, (ftnlen) 1))
	{
		i__1 = *k;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			if (tau[i__] == 0.)
			{

				/*              H(i)  =  I */

				i__2 = i__;
				for (j = 1; j <= i__2; ++j)
				{
					t[j + i__ * t_dim1] = 0.;
					/* L10: */
				}
			}
			else
			{

				/*              general case */

				vii = v[i__ + i__ * v_dim1];
				v[i__ + i__ * v_dim1] = 1.;
				if (lsame_ (storev, "C", (ftnlen) 1, (ftnlen) 1))
				{

					/*                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i) */

					i__2 = *n - i__ + 1;
					i__3 = i__ - 1;
					d__1 = -tau[i__];
					dgemv_ ("Transpose", &i__2, &i__3, &d__1, &v[i__ + v_dim1],
							  ldv, &v[i__ + i__ * v_dim1], &c__1, &c_b507, &t[i__ * t_dim1 + 1], &c__1, (ftnlen) 9);
				}
				else
				{

					/*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)' */

					i__2 = i__ - 1;
					i__3 = *n - i__ + 1;
					d__1 = -tau[i__];
					dgemv_ ("No transpose", &i__2, &i__3, &d__1, &v[i__ *
																					v_dim1 + 1], ldv, &v[i__ + i__ * v_dim1], ldv, &c_b507, &t[i__ * t_dim1 + 1], &c__1,
							  (ftnlen) 12);
				}
				v[i__ + i__ * v_dim1] = vii;

				/*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i) */

				i__2 = i__ - 1;
				dtrmv_ ("Upper", "No transpose", "Non-unit", &i__2, &t[t_offset], ldt, &t[i__ * t_dim1 + 1], &c__1, (ftnlen) 5, (ftnlen) 12, (ftnlen) 8);
				t[i__ + i__ * t_dim1] = tau[i__];
			}
			/* L20: */
		}
	}
	else
	{
		for (i__ = *k; i__ >= 1; --i__)
		{
			if (tau[i__] == 0.)
			{

				/*              H(i)  =  I */

				i__1 = *k;
				for (j = i__; j <= i__1; ++j)
				{
					t[j + i__ * t_dim1] = 0.;
					/* L30: */
				}
			}
			else
			{

				/*              general case */

				if (i__ < *k)
				{
					if (lsame_ (storev, "C", (ftnlen) 1, (ftnlen) 1))
					{
						vii = v[*n - *k + i__ + i__ * v_dim1];
						v[*n - *k + i__ + i__ * v_dim1] = 1.;

						/*                    T(i+1:k,i) := */
						/*                            - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i) */

						i__1 = *n - *k + i__;
						i__2 = *k - i__;
						d__1 = -tau[i__];
						dgemv_ ("Transpose", &i__1, &i__2, &d__1, &v[(i__ + 1)
																					* v_dim1 + 1], ldv, &v[i__ * v_dim1 + 1], &c__1, &c_b507, &t[i__ + 1 + i__ * t_dim1], &c__1,
								  (ftnlen) 9);
						v[*n - *k + i__ + i__ * v_dim1] = vii;
					}
					else
					{
						vii = v[i__ + (*n - *k + i__) * v_dim1];
						v[i__ + (*n - *k + i__) * v_dim1] = 1.;

						/*                    T(i+1:k,i) := */
						/*                            - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)' */

						i__1 = *k - i__;
						i__2 = *n - *k + i__;
						d__1 = -tau[i__];
						dgemv_ ("No transpose", &i__1, &i__2, &d__1, &v[i__ +
																						1 + v_dim1], ldv, &v[i__ + v_dim1], ldv, &c_b507, &t[i__ + 1 + i__ * t_dim1], &c__1,
								  (ftnlen) 12);
						v[i__ + (*n - *k + i__) * v_dim1] = vii;
					}

					/*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i) */

					i__1 = *k - i__;
					dtrmv_ ("Lower", "No transpose", "Non-unit", &i__1, &t[i__
																							 + 1 + (i__ + 1) * t_dim1], ldt, &t[i__ + 1 + i__ *
																																			t_dim1], &c__1, (ftnlen) 5, (ftnlen) 12, (ftnlen) 8);
				}
				t[i__ + i__ * t_dim1] = tau[i__];
			}
			/* L40: */
		}
	}
	return 0;

	/*     End of DLARFT */

}	/* dlarft_ */

static int
dlarfx_ (char *side, integer * m, integer * n, doublereal * v, doublereal * tau, doublereal * c__, integer * ldc, doublereal * work, ftnlen side_len)
{
	/* System generated locals */
	integer c_dim1, c_offset, i__1;
	doublereal d__1;

	/* Local variables */
	static integer j;
	static doublereal t1, t2, t3, t4, t5, t6, t7, t8, t9, v1, v2, v3, v4, v5, v6, v7, v8, v9, t10, v10, sum;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     February 29, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLARFX applies a real elementary reflector H to a real m by n */
	/*  matrix C, from either the left or the right. H is represented in the */
	/*  form */

	/*        H = I - tau * v * v' */

	/*  where tau is a real scalar and v is a real vector. */

	/*  If tau = 0, then H is taken to be the unit matrix */

	/*  This version uses inline code if H has order < 11. */

	/*  Arguments */
	/*  ========= */

	/*  SIDE    (input) CHARACTER*1 */
	/*          = 'L': form  H * C */
	/*          = 'R': form  C * H */

	/*  M       (input) INTEGER */
	/*          The number of rows of the matrix C. */

	/*  N       (input) INTEGER */
	/*          The number of columns of the matrix C. */

	/*  V       (input) DOUBLE PRECISION array, dimension (M) if SIDE = 'L' */
	/*                                     or (N) if SIDE = 'R' */
	/*          The vector v in the representation of H. */

	/*  TAU     (input) DOUBLE PRECISION */
	/*          The value tau in the representation of H. */

	/*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
	/*          On entry, the m by n matrix C. */
	/*          On exit, C is overwritten by the matrix H * C if SIDE = 'L', */
	/*          or C * H if SIDE = 'R'. */

	/*  LDC     (input) INTEGER */
	/*          The leading dimension of the array C. LDA >= (1,M). */

	/*  WORK    (workspace) DOUBLE PRECISION array, dimension */
	/*                      (N) if SIDE = 'L' */
	/*                      or (M) if SIDE = 'R' */
	/*          WORK is not referenced if H has order < 11. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/* Parameter adjustments */
	--v;
	c_dim1 = *ldc;
	c_offset = 1 + c_dim1;
	c__ -= c_offset;
	--work;

	/* Function Body */
	if (*tau == 0.)
	{
		return 0;
	}
	if (lsame_ (side, "L", (ftnlen) 1, (ftnlen) 1))
	{

		/*        Form  H * C, where H has order m. */

		switch (*m)
		{
		case 1:
			goto L10;
		case 2:
			goto L30;
		case 3:
			goto L50;
		case 4:
			goto L70;
		case 5:
			goto L90;
		case 6:
			goto L110;
		case 7:
			goto L130;
		case 8:
			goto L150;
		case 9:
			goto L170;
		case 10:
			goto L190;
		}

		/*        Code for general M */

		/*        w := C'*v */

		dgemv_ ("Transpose", m, n, &c_b348, &c__[c_offset], ldc, &v[1], &c__1, &c_b507, &work[1], &c__1, (ftnlen) 9);

		/*        C := C - tau * v * w' */

		d__1 = -(*tau);
		dger_ (m, n, &d__1, &v[1], &c__1, &work[1], &c__1, &c__[c_offset], ldc);
		goto L410;
	 L10:

		/*        Special code for 1 x 1 Householder */

		t1 = 1. - *tau * v[1] * v[1];
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			c__[j * c_dim1 + 1] = t1 * c__[j * c_dim1 + 1];
			/* L20: */
		}
		goto L410;
	 L30:

		/*        Special code for 2 x 2 Householder */

		v1 = v[1];
		t1 = *tau * v1;
		v2 = v[2];
		t2 = *tau * v2;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2];
			c__[j * c_dim1 + 1] -= sum * t1;
			c__[j * c_dim1 + 2] -= sum * t2;
			/* L40: */
		}
		goto L410;
	 L50:

		/*        Special code for 3 x 3 Householder */

		v1 = v[1];
		t1 = *tau * v1;
		v2 = v[2];
		t2 = *tau * v2;
		v3 = v[3];
		t3 = *tau * v3;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * c__[j * c_dim1 + 3];
			c__[j * c_dim1 + 1] -= sum * t1;
			c__[j * c_dim1 + 2] -= sum * t2;
			c__[j * c_dim1 + 3] -= sum * t3;
			/* L60: */
		}
		goto L410;
	 L70:

		/*        Special code for 4 x 4 Householder */

		v1 = v[1];
		t1 = *tau * v1;
		v2 = v[2];
		t2 = *tau * v2;
		v3 = v[3];
		t3 = *tau * v3;
		v4 = v[4];
		t4 = *tau * v4;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4];
			c__[j * c_dim1 + 1] -= sum * t1;
			c__[j * c_dim1 + 2] -= sum * t2;
			c__[j * c_dim1 + 3] -= sum * t3;
			c__[j * c_dim1 + 4] -= sum * t4;
			/* L80: */
		}
		goto L410;
	 L90:

		/*        Special code for 5 x 5 Householder */

		v1 = v[1];
		t1 = *tau * v1;
		v2 = v[2];
		t2 = *tau * v2;
		v3 = v[3];
		t3 = *tau * v3;
		v4 = v[4];
		t4 = *tau * v4;
		v5 = v[5];
		t5 = *tau * v5;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 * c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[j * c_dim1 + 5];
			c__[j * c_dim1 + 1] -= sum * t1;
			c__[j * c_dim1 + 2] -= sum * t2;
			c__[j * c_dim1 + 3] -= sum * t3;
			c__[j * c_dim1 + 4] -= sum * t4;
			c__[j * c_dim1 + 5] -= sum * t5;
			/* L100: */
		}
		goto L410;
	 L110:

		/*        Special code for 6 x 6 Householder */

		v1 = v[1];
		t1 = *tau * v1;
		v2 = v[2];
		t2 = *tau * v2;
		v3 = v[3];
		t3 = *tau * v3;
		v4 = v[4];
		t4 = *tau * v4;
		v5 = v[5];
		t5 = *tau * v5;
		v6 = v[6];
		t6 = *tau * v6;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 *
				c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6];
			c__[j * c_dim1 + 1] -= sum * t1;
			c__[j * c_dim1 + 2] -= sum * t2;
			c__[j * c_dim1 + 3] -= sum * t3;
			c__[j * c_dim1 + 4] -= sum * t4;
			c__[j * c_dim1 + 5] -= sum * t5;
			c__[j * c_dim1 + 6] -= sum * t6;
			/* L120: */
		}
		goto L410;
	 L130:

		/*        Special code for 7 x 7 Householder */

		v1 = v[1];
		t1 = *tau * v1;
		v2 = v[2];
		t2 = *tau * v2;
		v3 = v[3];
		t3 = *tau * v3;
		v4 = v[4];
		t4 = *tau * v4;
		v5 = v[5];
		t5 = *tau * v5;
		v6 = v[6];
		t6 = *tau * v6;
		v7 = v[7];
		t7 = *tau * v7;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 *
				c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6] + v7 * c__[j * c_dim1 + 7];
			c__[j * c_dim1 + 1] -= sum * t1;
			c__[j * c_dim1 + 2] -= sum * t2;
			c__[j * c_dim1 + 3] -= sum * t3;
			c__[j * c_dim1 + 4] -= sum * t4;
			c__[j * c_dim1 + 5] -= sum * t5;
			c__[j * c_dim1 + 6] -= sum * t6;
			c__[j * c_dim1 + 7] -= sum * t7;
			/* L140: */
		}
		goto L410;
	 L150:

		/*        Special code for 8 x 8 Householder */

		v1 = v[1];
		t1 = *tau * v1;
		v2 = v[2];
		t2 = *tau * v2;
		v3 = v[3];
		t3 = *tau * v3;
		v4 = v[4];
		t4 = *tau * v4;
		v5 = v[5];
		t5 = *tau * v5;
		v6 = v[6];
		t6 = *tau * v6;
		v7 = v[7];
		t7 = *tau * v7;
		v8 = v[8];
		t8 = *tau * v8;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 *
				c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6] + v7 * c__[j *
																																									 c_dim1 + 7] + v8 * c__[j * c_dim1 +
																																																	8];
			c__[j * c_dim1 + 1] -= sum * t1;
			c__[j * c_dim1 + 2] -= sum * t2;
			c__[j * c_dim1 + 3] -= sum * t3;
			c__[j * c_dim1 + 4] -= sum * t4;
			c__[j * c_dim1 + 5] -= sum * t5;
			c__[j * c_dim1 + 6] -= sum * t6;
			c__[j * c_dim1 + 7] -= sum * t7;
			c__[j * c_dim1 + 8] -= sum * t8;
			/* L160: */
		}
		goto L410;
	 L170:

		/*        Special code for 9 x 9 Householder */

		v1 = v[1];
		t1 = *tau * v1;
		v2 = v[2];
		t2 = *tau * v2;
		v3 = v[3];
		t3 = *tau * v3;
		v4 = v[4];
		t4 = *tau * v4;
		v5 = v[5];
		t5 = *tau * v5;
		v6 = v[6];
		t6 = *tau * v6;
		v7 = v[7];
		t7 = *tau * v7;
		v8 = v[8];
		t8 = *tau * v8;
		v9 = v[9];
		t9 = *tau * v9;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 *
				c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6] + v7 * c__[j *
																																									 c_dim1 + 7] + v8 * c__[j * c_dim1 +
																																																	8] +
				v9 * c__[j * c_dim1 + 9];
			c__[j * c_dim1 + 1] -= sum * t1;
			c__[j * c_dim1 + 2] -= sum * t2;
			c__[j * c_dim1 + 3] -= sum * t3;
			c__[j * c_dim1 + 4] -= sum * t4;
			c__[j * c_dim1 + 5] -= sum * t5;
			c__[j * c_dim1 + 6] -= sum * t6;
			c__[j * c_dim1 + 7] -= sum * t7;
			c__[j * c_dim1 + 8] -= sum * t8;
			c__[j * c_dim1 + 9] -= sum * t9;
			/* L180: */
		}
		goto L410;
	 L190:

		/*        Special code for 10 x 10 Householder */

		v1 = v[1];
		t1 = *tau * v1;
		v2 = v[2];
		t2 = *tau * v2;
		v3 = v[3];
		t3 = *tau * v3;
		v4 = v[4];
		t4 = *tau * v4;
		v5 = v[5];
		t5 = *tau * v5;
		v6 = v[6];
		t6 = *tau * v6;
		v7 = v[7];
		t7 = *tau * v7;
		v8 = v[8];
		t8 = *tau * v8;
		v9 = v[9];
		t9 = *tau * v9;
		v10 = v[10];
		t10 = *tau * v10;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			sum = v1 * c__[j * c_dim1 + 1] + v2 * c__[j * c_dim1 + 2] + v3 *
				c__[j * c_dim1 + 3] + v4 * c__[j * c_dim1 + 4] + v5 * c__[j * c_dim1 + 5] + v6 * c__[j * c_dim1 + 6] + v7 * c__[j *
																																									 c_dim1 + 7] + v8 * c__[j * c_dim1 +
																																																	8] +
				v9 * c__[j * c_dim1 + 9] + v10 * c__[j * c_dim1 + 10];
			c__[j * c_dim1 + 1] -= sum * t1;
			c__[j * c_dim1 + 2] -= sum * t2;
			c__[j * c_dim1 + 3] -= sum * t3;
			c__[j * c_dim1 + 4] -= sum * t4;
			c__[j * c_dim1 + 5] -= sum * t5;
			c__[j * c_dim1 + 6] -= sum * t6;
			c__[j * c_dim1 + 7] -= sum * t7;
			c__[j * c_dim1 + 8] -= sum * t8;
			c__[j * c_dim1 + 9] -= sum * t9;
			c__[j * c_dim1 + 10] -= sum * t10;
			/* L200: */
		}
		goto L410;
	}
	else
	{

		/*        Form  C * H, where H has order n. */

		switch (*n)
		{
		case 1:
			goto L210;
		case 2:
			goto L230;
		case 3:
			goto L250;
		case 4:
			goto L270;
		case 5:
			goto L290;
		case 6:
			goto L310;
		case 7:
			goto L330;
		case 8:
			goto L350;
		case 9:
			goto L370;
		case 10:
			goto L390;
		}

		/*        Code for general N */

		/*        w := C * v */

		dgemv_ ("No transpose", m, n, &c_b348, &c__[c_offset], ldc, &v[1], &c__1, &c_b507, &work[1], &c__1, (ftnlen) 12);

		/*        C := C - tau * w * v' */

		d__1 = -(*tau);
		dger_ (m, n, &d__1, &work[1], &c__1, &v[1], &c__1, &c__[c_offset], ldc);
		goto L410;
	 L210:

		/*        Special code for 1 x 1 Householder */

		t1 = 1. - *tau * v[1] * v[1];
		i__1 = *m;
		for (j = 1; j <= i__1; ++j)
		{
			c__[j + c_dim1] = t1 * c__[j + c_dim1];
			/* L220: */
		}
		goto L410;
	 L230:

		/*        Special code for 2 x 2 Householder */

		v1 = v[1];
		t1 = *tau * v1;
		v2 = v[2];
		t2 = *tau * v2;
		i__1 = *m;
		for (j = 1; j <= i__1; ++j)
		{
			sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)];
			c__[j + c_dim1] -= sum * t1;
			c__[j + (c_dim1 << 1)] -= sum * t2;
			/* L240: */
		}
		goto L410;
	 L250:

		/*        Special code for 3 x 3 Householder */

		v1 = v[1];
		t1 = *tau * v1;
		v2 = v[2];
		t2 = *tau * v2;
		v3 = v[3];
		t3 = *tau * v3;
		i__1 = *m;
		for (j = 1; j <= i__1; ++j)
		{
			sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * c__[j + c_dim1 * 3];
			c__[j + c_dim1] -= sum * t1;
			c__[j + (c_dim1 << 1)] -= sum * t2;
			c__[j + c_dim1 * 3] -= sum * t3;
			/* L260: */
		}
		goto L410;
	 L270:

		/*        Special code for 4 x 4 Householder */

		v1 = v[1];
		t1 = *tau * v1;
		v2 = v[2];
		t2 = *tau * v2;
		v3 = v[3];
		t3 = *tau * v3;
		v4 = v[4];
		t4 = *tau * v4;
		i__1 = *m;
		for (j = 1; j <= i__1; ++j)
		{
			sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)];
			c__[j + c_dim1] -= sum * t1;
			c__[j + (c_dim1 << 1)] -= sum * t2;
			c__[j + c_dim1 * 3] -= sum * t3;
			c__[j + (c_dim1 << 2)] -= sum * t4;
			/* L280: */
		}
		goto L410;
	 L290:

		/*        Special code for 5 x 5 Householder */

		v1 = v[1];
		t1 = *tau * v1;
		v2 = v[2];
		t2 = *tau * v2;
		v3 = v[3];
		t3 = *tau * v3;
		v4 = v[4];
		t4 = *tau * v4;
		v5 = v[5];
		t5 = *tau * v5;
		i__1 = *m;
		for (j = 1; j <= i__1; ++j)
		{
			sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 * c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * c__[j + c_dim1 * 5];
			c__[j + c_dim1] -= sum * t1;
			c__[j + (c_dim1 << 1)] -= sum * t2;
			c__[j + c_dim1 * 3] -= sum * t3;
			c__[j + (c_dim1 << 2)] -= sum * t4;
			c__[j + c_dim1 * 5] -= sum * t5;
			/* L300: */
		}
		goto L410;
	 L310:

		/*        Special code for 6 x 6 Householder */

		v1 = v[1];
		t1 = *tau * v1;
		v2 = v[2];
		t2 = *tau * v2;
		v3 = v[3];
		t3 = *tau * v3;
		v4 = v[4];
		t4 = *tau * v4;
		v5 = v[5];
		t5 = *tau * v5;
		v6 = v[6];
		t6 = *tau * v6;
		i__1 = *m;
		for (j = 1; j <= i__1; ++j)
		{
			sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 *
				c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6];
			c__[j + c_dim1] -= sum * t1;
			c__[j + (c_dim1 << 1)] -= sum * t2;
			c__[j + c_dim1 * 3] -= sum * t3;
			c__[j + (c_dim1 << 2)] -= sum * t4;
			c__[j + c_dim1 * 5] -= sum * t5;
			c__[j + c_dim1 * 6] -= sum * t6;
			/* L320: */
		}
		goto L410;
	 L330:

		/*        Special code for 7 x 7 Householder */

		v1 = v[1];
		t1 = *tau * v1;
		v2 = v[2];
		t2 = *tau * v2;
		v3 = v[3];
		t3 = *tau * v3;
		v4 = v[4];
		t4 = *tau * v4;
		v5 = v[5];
		t5 = *tau * v5;
		v6 = v[6];
		t6 = *tau * v6;
		v7 = v[7];
		t7 = *tau * v7;
		i__1 = *m;
		for (j = 1; j <= i__1; ++j)
		{
			sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 *
				c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 * c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6] + v7 * c__[j + c_dim1 * 7];
			c__[j + c_dim1] -= sum * t1;
			c__[j + (c_dim1 << 1)] -= sum * t2;
			c__[j + c_dim1 * 3] -= sum * t3;
			c__[j + (c_dim1 << 2)] -= sum * t4;
			c__[j + c_dim1 * 5] -= sum * t5;
			c__[j + c_dim1 * 6] -= sum * t6;
			c__[j + c_dim1 * 7] -= sum * t7;
			/* L340: */
		}
		goto L410;
	 L350:

		/*        Special code for 8 x 8 Householder */

		v1 = v[1];
		t1 = *tau * v1;
		v2 = v[2];
		t2 = *tau * v2;
		v3 = v[3];
		t3 = *tau * v3;
		v4 = v[4];
		t4 = *tau * v4;
		v5 = v[5];
		t5 = *tau * v5;
		v6 = v[6];
		t6 = *tau * v6;
		v7 = v[7];
		t7 = *tau * v7;
		v8 = v[8];
		t8 = *tau * v8;
		i__1 = *m;
		for (j = 1; j <= i__1; ++j)
		{
			sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 *
				c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 *
				c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6] + v7 * c__[j + c_dim1 * 7] + v8 * c__[j + (c_dim1 << 3)];
			c__[j + c_dim1] -= sum * t1;
			c__[j + (c_dim1 << 1)] -= sum * t2;
			c__[j + c_dim1 * 3] -= sum * t3;
			c__[j + (c_dim1 << 2)] -= sum * t4;
			c__[j + c_dim1 * 5] -= sum * t5;
			c__[j + c_dim1 * 6] -= sum * t6;
			c__[j + c_dim1 * 7] -= sum * t7;
			c__[j + (c_dim1 << 3)] -= sum * t8;
			/* L360: */
		}
		goto L410;
	 L370:

		/*        Special code for 9 x 9 Householder */

		v1 = v[1];
		t1 = *tau * v1;
		v2 = v[2];
		t2 = *tau * v2;
		v3 = v[3];
		t3 = *tau * v3;
		v4 = v[4];
		t4 = *tau * v4;
		v5 = v[5];
		t5 = *tau * v5;
		v6 = v[6];
		t6 = *tau * v6;
		v7 = v[7];
		t7 = *tau * v7;
		v8 = v[8];
		t8 = *tau * v8;
		v9 = v[9];
		t9 = *tau * v9;
		i__1 = *m;
		for (j = 1; j <= i__1; ++j)
		{
			sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 *
				c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 *
				c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6] + v7 * c__[j + c_dim1 * 7] + v8 * c__[j + (c_dim1 << 3)] + v9 * c__[j + c_dim1 * 9];
			c__[j + c_dim1] -= sum * t1;
			c__[j + (c_dim1 << 1)] -= sum * t2;
			c__[j + c_dim1 * 3] -= sum * t3;
			c__[j + (c_dim1 << 2)] -= sum * t4;
			c__[j + c_dim1 * 5] -= sum * t5;
			c__[j + c_dim1 * 6] -= sum * t6;
			c__[j + c_dim1 * 7] -= sum * t7;
			c__[j + (c_dim1 << 3)] -= sum * t8;
			c__[j + c_dim1 * 9] -= sum * t9;
			/* L380: */
		}
		goto L410;
	 L390:

		/*        Special code for 10 x 10 Householder */

		v1 = v[1];
		t1 = *tau * v1;
		v2 = v[2];
		t2 = *tau * v2;
		v3 = v[3];
		t3 = *tau * v3;
		v4 = v[4];
		t4 = *tau * v4;
		v5 = v[5];
		t5 = *tau * v5;
		v6 = v[6];
		t6 = *tau * v6;
		v7 = v[7];
		t7 = *tau * v7;
		v8 = v[8];
		t8 = *tau * v8;
		v9 = v[9];
		t9 = *tau * v9;
		v10 = v[10];
		t10 = *tau * v10;
		i__1 = *m;
		for (j = 1; j <= i__1; ++j)
		{
			sum = v1 * c__[j + c_dim1] + v2 * c__[j + (c_dim1 << 1)] + v3 *
				c__[j + c_dim1 * 3] + v4 * c__[j + (c_dim1 << 2)] + v5 *
				c__[j + c_dim1 * 5] + v6 * c__[j + c_dim1 * 6] + v7 * c__[j + c_dim1 * 7] + v8 * c__[j + (c_dim1 << 3)] + v9 * c__[j + c_dim1 * 9] + v10 * c__[j +
																																																			  c_dim1
																																																			  *
																																																			  10];
			c__[j + c_dim1] -= sum * t1;
			c__[j + (c_dim1 << 1)] -= sum * t2;
			c__[j + c_dim1 * 3] -= sum * t3;
			c__[j + (c_dim1 << 2)] -= sum * t4;
			c__[j + c_dim1 * 5] -= sum * t5;
			c__[j + c_dim1 * 6] -= sum * t6;
			c__[j + c_dim1 * 7] -= sum * t7;
			c__[j + (c_dim1 << 3)] -= sum * t8;
			c__[j + c_dim1 * 9] -= sum * t9;
			c__[j + c_dim1 * 10] -= sum * t10;
			/* L400: */
		}
		goto L410;
	}
 L410:
	return 0;

	/*     End of DLARFX */

}	/* dlarfx_ */

static int
dlarnv_ (integer * idist, integer * iseed, integer * n, doublereal * x)
{
	/* System generated locals */
	integer i__1, i__2, i__3;

	/* Builtin functions */
	double log (doublereal), sqrt (doublereal), cos (doublereal);

	/* Local variables */
	static integer i__;
	static doublereal u[128];
	static integer il, iv, il2;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     September 30, 1994 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLARNV returns a vector of n random real numbers from a uniform or */
	/*  normal distribution. */

	/*  Arguments */
	/*  ========= */

	/*  IDIST   (input) INTEGER */
	/*          Specifies the distribution of the random numbers: */
	/*          = 1:  uniform (0,1) */
	/*          = 2:  uniform (-1,1) */
	/*          = 3:  normal (0,1) */

	/*  ISEED   (input/output) INTEGER array, dimension (4) */
	/*          On entry, the seed of the random number generator; the array */
	/*          elements must be between 0 and 4095, and ISEED(4) must be */
	/*          odd. */
	/*          On exit, the seed is updated. */

	/*  N       (input) INTEGER */
	/*          The number of random numbers to be generated. */

	/*  X       (output) DOUBLE PRECISION array, dimension (N) */
	/*          The generated random numbers. */

	/*  Further Details */
	/*  =============== */

	/*  This routine calls the auxiliary routine DLARUV to generate random */
	/*  real numbers from a uniform (0,1) distribution, in batches of up to */
	/*  128 using vectorisable code. The Box-Muller method is used to */
	/*  transform numbers from a uniform to a normal distribution. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Local Arrays .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/* Parameter adjustments */
	--x;
	--iseed;

	/* Function Body */
	i__1 = *n;
	for (iv = 1; iv <= i__1; iv += 64)
	{
		/* Computing MIN */
		i__2 = 64, i__3 = *n - iv + 1;
		il = min (i__2, i__3);
		if (*idist == 3)
		{
			il2 = il << 1;
		}
		else
		{
			il2 = il;
		}

		/*        Call DLARUV to generate IL2 numbers from a uniform (0,1) */
		/*        distribution (IL2 <= LV) */

		dlaruv_ (&iseed[1], &il2, u);

		if (*idist == 1)
		{

			/*           Copy generated numbers */

			i__2 = il;
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				x[iv + i__ - 1] = u[i__ - 1];
				/* L10: */
			}
		}
		else if (*idist == 2)
		{

			/*           Convert generated numbers to uniform (-1,1) distribution */

			i__2 = il;
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				x[iv + i__ - 1] = u[i__ - 1] * 2. - 1.;
				/* L20: */
			}
		}
		else if (*idist == 3)
		{

			/*           Convert generated numbers to normal (0,1) distribution */

			i__2 = il;
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				x[iv + i__ - 1] = sqrt (log (u[(i__ << 1) - 2]) * -2.) * cos (u[(i__ << 1) - 1] * 6.2831853071795864769252867663);
				/* L30: */
			}
		}
		/* L40: */
	}
	return 0;

	/*     End of DLARNV */

}	/* dlarnv_ */

static int
dlartg_ (doublereal * f, doublereal * g, doublereal * cs, doublereal * sn, doublereal * r__)
{
	/* Initialized data */

	static logical first = TRUE_;

	/* System generated locals */
	integer i__1;
	doublereal d__1, d__2;

	/* Builtin functions */
	double log (doublereal), pow_di (doublereal *, integer *), sqrt (doublereal);

	/* Local variables */
	static integer i__;
	static doublereal f1, g1, eps, scale;
	static integer count;
	static doublereal safmn2, safmx2;
	static doublereal safmin;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     September 30, 1994 */

	/*     .. Scalar Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLARTG generate a plane rotation so that */

	/*     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1. */
	/*     [ -SN  CS  ]     [ G ]     [ 0 ] */

	/*  This is a slower, more accurate version of the BLAS1 routine DROTG, */
	/*  with the following other differences: */
	/*     F and G are unchanged on return. */
	/*     If G=0, then CS=1 and SN=0. */
	/*     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any */
	/*        floating point operations (saves work in DBDSQR when */
	/*        there are zeros on the diagonal). */

	/*  If F exceeds G in magnitude, CS will be positive. */

	/*  Arguments */
	/*  ========= */

	/*  F       (input) DOUBLE PRECISION */
	/*          The first component of vector to be rotated. */

	/*  G       (input) DOUBLE PRECISION */
	/*          The second component of vector to be rotated. */

	/*  CS      (output) DOUBLE PRECISION */
	/*          The cosine of the rotation. */

	/*  SN      (output) DOUBLE PRECISION */
	/*          The sine of the rotation. */

	/*  R       (output) DOUBLE PRECISION */
	/*          The nonzero component of the rotated vector. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Save statement .. */
	/*     .. */
	/*     .. Data statements .. */
	/*     .. */
	/*     .. Executable Statements .. */

	if (first)
	{
		first = FALSE_;
		safmin = dlamch_ ("S", (ftnlen) 1);
		eps = dlamch_ ("E", (ftnlen) 1);
		d__1 = dlamch_ ("B", (ftnlen) 1);
		i__1 = (integer) (log (safmin / eps) / log (dlamch_ ("B", (ftnlen) 1)) / 2.);
		safmn2 = pow_di (&d__1, &i__1);
		safmx2 = 1. / safmn2;
	}
	if (*g == 0.)
	{
		*cs = 1.;
		*sn = 0.;
		*r__ = *f;
	}
	else if (*f == 0.)
	{
		*cs = 0.;
		*sn = 1.;
		*r__ = *g;
	}
	else
	{
		f1 = *f;
		g1 = *g;
		/* Computing MAX */
		d__1 = abs (f1), d__2 = abs (g1);
		scale = max (d__1, d__2);
		if (scale >= safmx2)
		{
			count = 0;
		 L10:
			++count;
			f1 *= safmn2;
			g1 *= safmn2;
			/* Computing MAX */
			d__1 = abs (f1), d__2 = abs (g1);
			scale = max (d__1, d__2);
			if (scale >= safmx2)
			{
				goto L10;
			}
			/* Computing 2nd power */
			d__1 = f1;
			/* Computing 2nd power */
			d__2 = g1;
			*r__ = sqrt (d__1 * d__1 + d__2 * d__2);
			*cs = f1 / *r__;
			*sn = g1 / *r__;
			i__1 = count;
			for (i__ = 1; i__ <= i__1; ++i__)
			{
				*r__ *= safmx2;
				/* L20: */
			}
		}
		else if (scale <= safmn2)
		{
			count = 0;
		 L30:
			++count;
			f1 *= safmx2;
			g1 *= safmx2;
			/* Computing MAX */
			d__1 = abs (f1), d__2 = abs (g1);
			scale = max (d__1, d__2);
			if (scale <= safmn2)
			{
				goto L30;
			}
			/* Computing 2nd power */
			d__1 = f1;
			/* Computing 2nd power */
			d__2 = g1;
			*r__ = sqrt (d__1 * d__1 + d__2 * d__2);
			*cs = f1 / *r__;
			*sn = g1 / *r__;
			i__1 = count;
			for (i__ = 1; i__ <= i__1; ++i__)
			{
				*r__ *= safmn2;
				/* L40: */
			}
		}
		else
		{
			/* Computing 2nd power */
			d__1 = f1;
			/* Computing 2nd power */
			d__2 = g1;
			*r__ = sqrt (d__1 * d__1 + d__2 * d__2);
			*cs = f1 / *r__;
			*sn = g1 / *r__;
		}
		if (abs (*f) > abs (*g) && *cs < 0.)
		{
			*cs = -(*cs);
			*sn = -(*sn);
			*r__ = -(*r__);
		}
	}
	return 0;

	/*     End of DLARTG */

}	/* dlartg_ */

static int
dlaruv_ (integer * iseed, integer * n, doublereal * x)
{
	/* Initialized data */

	static integer mm[512] /* was [128][4] */  = { 494, 2637, 255, 2008, 1253,
		3344, 4084, 1739, 3143, 3468, 688, 1657, 1238, 3166, 1292, 3422, 1270, 2016,
		154, 2862, 697, 1706, 491, 931, 1444, 444, 3577, 3944, 2184, 1661, 3482, 657,
		3023, 3618, 1267, 1828, 164, 3798, 3087, 2400, 2870, 3876, 1905, 1593, 1797,
		1234, 3460, 328, 2861, 1950, 617, 2070, 3331, 769, 1558, 2412, 2800, 189, 287,
		2045, 1227, 2838, 209, 2770, 3654, 3993, 192, 2253, 3491, 2889, 2857, 2094,
		1818, 688, 1407, 634, 3231, 815, 3524, 1914, 516, 164, 303, 2144, 3480, 119,
		3357, 837, 2826, 2332, 2089, 3780, 1700, 3712, 150, 2000, 3375, 1621, 3090,
		3765, 1149, 3146, 33, 3082, 2741, 359, 3316, 1749, 185, 2784, 2202, 2199, 1364,
		1244, 2020, 3160, 2785, 2772, 1217, 1822, 1245, 2252, 3904, 2774, 997, 2573,
		1148, 545, 322, 789, 1440, 752, 2859, 123, 1848, 643, 2405, 2638, 2344, 46,
		3814, 913, 3649, 339, 3808, 822, 2832, 3078, 3633, 2970, 637, 2249, 2081, 4019,
		1478, 242, 481, 2075, 4058, 622, 3376, 812, 234, 641, 4005, 1122, 3135, 2640,
		2302, 40, 1832, 2247, 2034, 2637, 1287, 1691, 496, 1597, 2394, 2584, 1843, 336,
		1472, 2407, 433, 2096, 1761, 2810, 566, 442, 41, 1238, 1086, 603, 840, 3168,
		1499, 1084, 3438, 2408, 1589, 2391, 288, 26, 512, 1456, 171, 1677, 2657, 2270,
		2587, 2961, 1970, 1817, 676, 1410, 3723, 2803, 3185, 184, 663, 499, 3784, 1631,
		1925, 3912, 1398, 1349, 1441, 2224, 2411, 1907, 3192, 2786, 382, 37, 759, 2948,
		1862, 3802, 2423, 2051, 2295, 1332, 1832, 2405, 3638, 3661, 327, 3660, 716,
		1842, 3987, 1368, 1848, 2366, 2508, 3754, 1766, 3572, 2893, 307, 1297, 3966,
		758, 2598, 3406, 2922, 1038, 2934, 2091, 2451, 1580, 1958, 2055, 1507, 1078,
		3273, 17, 854, 2916, 3971, 2889, 3831, 2621, 1541, 893, 736, 3992, 787, 2125,
		2364, 2460, 257, 1574, 3912, 1216, 3248, 3401, 2124, 2762, 149, 2245, 166, 466,
		4018, 1399, 190, 2879, 153, 2320, 18, 712, 2159, 2318, 2091, 3443, 1510, 449,
		1956, 2201, 3137, 3399, 1321, 2271, 3667, 2703, 629, 2365, 2431, 1113, 3922,
		2554, 184, 2099, 3228, 4012, 1921, 3452, 3901, 572, 3309, 3171, 817, 3039,
		1696, 1256, 3715, 2077, 3019, 1497, 1101, 717, 51, 981, 1978, 1813, 3881, 76,
		3846, 3694, 1682, 124, 1660, 3997, 479, 1141, 886, 3514, 1301, 3604, 1888,
		1836, 1990, 2058, 692, 1194, 20, 3285, 2046, 2107, 3508, 3525, 3801, 2549,
		1145, 2253, 305, 3301, 1065, 3133, 2913, 3285, 1241, 1197, 3729, 2501, 1673,
		541, 2753, 949, 2361, 1165, 4081, 2725, 3305, 3069, 3617, 3733, 409, 2157,
		1361, 3973, 1865, 2525, 1409, 3445, 3577, 77, 3761, 2149, 1449, 3005, 225, 85,
		3673, 3117, 3089, 1349, 2057, 413, 65, 1845, 697, 3085, 3441, 1573, 3689, 2941,
		929, 533, 2841, 4077, 721, 2821, 2249, 2397, 2817, 245, 1913, 1997, 3121, 997,
		1833, 2877, 1633, 981, 2009, 941, 2449, 197, 2441, 285, 1473, 2741, 3129, 909,
		2801, 421, 4073, 2813, 2337, 1429, 1177, 1901, 81, 1669, 2633, 2269, 129, 1141,
		249, 3917, 2481, 3941, 2217, 2749, 3041, 1877, 345, 2861, 1809, 3141, 2825,
		157, 2881, 3637, 1465, 2829, 2161, 3365, 361, 2685, 3745, 2325, 3609, 3821,
		3537, 517, 3017, 2141, 1537
	};

	/* System generated locals */
	integer i__1;

	/* Local variables */
	static integer i__, i1, i2, i3, i4, it1, it2, it3, it4;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     October 31, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLARUV returns a vector of n random real numbers from a uniform (0,1) */
	/*  distribution (n <= 128). */

	/*  This is an auxiliary routine called by DLARNV and ZLARNV. */

	/*  Arguments */
	/*  ========= */

	/*  ISEED   (input/output) INTEGER array, dimension (4) */
	/*          On entry, the seed of the random number generator; the array */
	/*          elements must be between 0 and 4095, and ISEED(4) must be */
	/*          odd. */
	/*          On exit, the seed is updated. */

	/*  N       (input) INTEGER */
	/*          The number of random numbers to be generated. N <= 128. */

	/*  X       (output) DOUBLE PRECISION array, dimension (N) */
	/*          The generated random numbers. */

	/*  Further Details */
	/*  =============== */

	/*  This routine uses a multiplicative congruential method with modulus */
	/*  2**48 and multiplier 33952834046453 (see G.S.Fishman, */
	/*  'Multiplicative congruential random number generators with modulus */
	/*  2**b: an exhaustive analysis for b = 32 and a partial analysis for */
	/*  b = 48', Math. Comp. 189, pp 331-344, 1990). */

	/*  48-bit integers are stored in 4 integer array elements with 12 bits */
	/*  per element. Hence the routine is portable across machines with */
	/*  integers of 32 bits or more. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Local Arrays .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Data statements .. */
	/* Parameter adjustments */
	--iseed;
	--x;

	/* Function Body */
	/*     .. */
	/*     .. Executable Statements .. */

	i1 = iseed[1];
	i2 = iseed[2];
	i3 = iseed[3];
	i4 = iseed[4];

	i__1 = min (*n, 128);
	for (i__ = 1; i__ <= i__1; ++i__)
	{

		/*        Multiply the seed by i-th power of the multiplier modulo 2**48 */

		it4 = i4 * mm[i__ + 383];
		it3 = it4 / 4096;
		it4 -= it3 << 12;
		it3 = it3 + i3 * mm[i__ + 383] + i4 * mm[i__ + 255];
		it2 = it3 / 4096;
		it3 -= it2 << 12;
		it2 = it2 + i2 * mm[i__ + 383] + i3 * mm[i__ + 255] + i4 * mm[i__ + 127];
		it1 = it2 / 4096;
		it2 -= it1 << 12;
		it1 = it1 + i1 * mm[i__ + 383] + i2 * mm[i__ + 255] + i3 * mm[i__ + 127] + i4 * mm[i__ - 1];
		it1 %= 4096;

		/*        Convert 48-bit integer to a real number in the interval (0,1) */

		x[i__] = ((doublereal) it1 + ((doublereal) it2 + ((doublereal) it3 + (doublereal) it4 * 2.44140625e-4) * 2.44140625e-4) * 2.44140625e-4) * 2.44140625e-4;
		/* L10: */
	}

	/*     Return final value of seed */

	iseed[1] = it1;
	iseed[2] = it2;
	iseed[3] = it3;
	iseed[4] = it4;
	return 0;

	/*     End of DLARUV */

}	/* dlaruv_ */

static int
dlascl_ (char *type__, integer * kl, integer * ku,
			doublereal * cfrom, doublereal * cto, integer * m, integer * n, doublereal * a, integer * lda, integer * info, ftnlen type_len)
{
	/* System generated locals */
	integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;

	/* Local variables */
	static integer i__, j, k1, k2, k3, k4;
	static doublereal mul, cto1;
	static logical done;
	static doublereal ctoc;
	static integer itype;
	static doublereal cfrom1;
	static doublereal cfromc;
	static doublereal bignum, smlnum;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     February 29, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLASCL multiplies the M by N real matrix A by the real scalar */
	/*  CTO/CFROM.  This is done without over/underflow as long as the final */
	/*  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that */
	/*  A may be full, upper triangular, lower triangular, upper Hessenberg, */
	/*  or banded. */

	/*  Arguments */
	/*  ========= */

	/*  TYPE    (input) CHARACTER*1 */
	/*          TYPE indices the storage type of the input matrix. */
	/*          = 'G':  A is a full matrix. */
	/*          = 'L':  A is a lower triangular matrix. */
	/*          = 'U':  A is an upper triangular matrix. */
	/*          = 'H':  A is an upper Hessenberg matrix. */
	/*          = 'B':  A is a symmetric band matrix with lower bandwidth KL */
	/*                  and upper bandwidth KU and with the only the lower */
	/*                  half stored. */
	/*          = 'Q':  A is a symmetric band matrix with lower bandwidth KL */
	/*                  and upper bandwidth KU and with the only the upper */
	/*                  half stored. */
	/*          = 'Z':  A is a band matrix with lower bandwidth KL and upper */
	/*                  bandwidth KU. */

	/*  KL      (input) INTEGER */
	/*          The lower bandwidth of A.  Referenced only if TYPE = 'B', */
	/*          'Q' or 'Z'. */

	/*  KU      (input) INTEGER */
	/*          The upper bandwidth of A.  Referenced only if TYPE = 'B', */
	/*          'Q' or 'Z'. */

	/*  CFROM   (input) DOUBLE PRECISION */
	/*  CTO     (input) DOUBLE PRECISION */
	/*          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed */
	/*          without over/underflow if the final result CTO*A(I,J)/CFROM */
	/*          can be represented without over/underflow.  CFROM must be */
	/*          nonzero. */

	/*  M       (input) INTEGER */
	/*          The number of rows of the matrix A.  M >= 0. */

	/*  N       (input) INTEGER */
	/*          The number of columns of the matrix A.  N >= 0. */

	/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,M) */
	/*          The matrix to be multiplied by CTO/CFROM.  See TYPE for the */
	/*          storage type. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(1,M). */

	/*  INFO    (output) INTEGER */
	/*          0  - successful exit */
	/*          <0 - if INFO = -i, the i-th argument had an illegal value. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Test the input arguments */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;

	/* Function Body */
	*info = 0;

	if (lsame_ (type__, "G", (ftnlen) 1, (ftnlen) 1))
	{
		itype = 0;
	}
	else if (lsame_ (type__, "L", (ftnlen) 1, (ftnlen) 1))
	{
		itype = 1;
	}
	else if (lsame_ (type__, "U", (ftnlen) 1, (ftnlen) 1))
	{
		itype = 2;
	}
	else if (lsame_ (type__, "H", (ftnlen) 1, (ftnlen) 1))
	{
		itype = 3;
	}
	else if (lsame_ (type__, "B", (ftnlen) 1, (ftnlen) 1))
	{
		itype = 4;
	}
	else if (lsame_ (type__, "Q", (ftnlen) 1, (ftnlen) 1))
	{
		itype = 5;
	}
	else if (lsame_ (type__, "Z", (ftnlen) 1, (ftnlen) 1))
	{
		itype = 6;
	}
	else
	{
		itype = -1;
	}

	if (itype == -1)
	{
		*info = -1;
	}
	else if (*cfrom == 0.)
	{
		*info = -4;
	}
	else if (*m < 0)
	{
		*info = -6;
	}
	else if (*n < 0 || itype == 4 && *n != *m || itype == 5 && *n != *m)
	{
		*info = -7;
	}
	else if (itype <= 3 && *lda < max (1, *m))
	{
		*info = -9;
	}
	else if (itype >= 4)
	{
		/* Computing MAX */
		i__1 = *m - 1;
		if (*kl < 0 || *kl > max (i__1, 0))
		{
			*info = -2;
		}
		else	/* if(complicated condition) */
		{
			/* Computing MAX */
			i__1 = *n - 1;
			if (*ku < 0 || *ku > max (i__1, 0) || (itype == 4 || itype == 5) && *kl != *ku)
			{
				*info = -3;
			}
			else if (itype == 4 && *lda < *kl + 1 || itype == 5 && *lda < *ku + 1 || itype == 6 && *lda < (*kl << 1) + *ku + 1)
			{
				*info = -9;
			}
		}
	}

	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DLASCL", &i__1, (ftnlen) 6);
		return 0;
	}

	/*     Quick return if possible */

	if (*n == 0 || *m == 0)
	{
		return 0;
	}

	/*     Get machine parameters */

	smlnum = dlamch_ ("S", (ftnlen) 1);
	bignum = 1. / smlnum;

	cfromc = *cfrom;
	ctoc = *cto;

 L10:
	cfrom1 = cfromc * smlnum;
	cto1 = ctoc / bignum;
	if (abs (cfrom1) > abs (ctoc) && ctoc != 0.)
	{
		mul = smlnum;
		done = FALSE_;
		cfromc = cfrom1;
	}
	else if (abs (cto1) > abs (cfromc))
	{
		mul = bignum;
		done = FALSE_;
		ctoc = cto1;
	}
	else
	{
		mul = ctoc / cfromc;
		done = TRUE_;
	}

	if (itype == 0)
	{

		/*        Full matrix */

		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				a[i__ + j * a_dim1] *= mul;
				/* L20: */
			}
			/* L30: */
		}

	}
	else if (itype == 1)
	{

		/*        Lower triangular matrix */

		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			i__2 = *m;
			for (i__ = j; i__ <= i__2; ++i__)
			{
				a[i__ + j * a_dim1] *= mul;
				/* L40: */
			}
			/* L50: */
		}

	}
	else if (itype == 2)
	{

		/*        Upper triangular matrix */

		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			i__2 = min (j, *m);
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				a[i__ + j * a_dim1] *= mul;
				/* L60: */
			}
			/* L70: */
		}

	}
	else if (itype == 3)
	{

		/*        Upper Hessenberg matrix */

		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			/* Computing MIN */
			i__3 = j + 1;
			i__2 = min (i__3, *m);
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				a[i__ + j * a_dim1] *= mul;
				/* L80: */
			}
			/* L90: */
		}

	}
	else if (itype == 4)
	{

		/*        Lower half of a symmetric band matrix */

		k3 = *kl + 1;
		k4 = *n + 1;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			/* Computing MIN */
			i__3 = k3, i__4 = k4 - j;
			i__2 = min (i__3, i__4);
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				a[i__ + j * a_dim1] *= mul;
				/* L100: */
			}
			/* L110: */
		}

	}
	else if (itype == 5)
	{

		/*        Upper half of a symmetric band matrix */

		k1 = *ku + 2;
		k3 = *ku + 1;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			/* Computing MAX */
			i__2 = k1 - j;
			i__3 = k3;
			for (i__ = max (i__2, 1); i__ <= i__3; ++i__)
			{
				a[i__ + j * a_dim1] *= mul;
				/* L120: */
			}
			/* L130: */
		}

	}
	else if (itype == 6)
	{

		/*        Band matrix */

		k1 = *kl + *ku + 2;
		k2 = *kl + 1;
		k3 = (*kl << 1) + *ku + 1;
		k4 = *kl + *ku + 1 + *m;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			/* Computing MAX */
			i__3 = k1 - j;
			/* Computing MIN */
			i__4 = k3, i__5 = k4 - j;
			i__2 = min (i__4, i__5);
			for (i__ = max (i__3, k2); i__ <= i__2; ++i__)
			{
				a[i__ + j * a_dim1] *= mul;
				/* L140: */
			}
			/* L150: */
		}

	}

	if (!done)
	{
		goto L10;
	}

	return 0;

	/*     End of DLASCL */

}	/* dlascl_ */

static int
dlaset_ (char *uplo, integer * m, integer * n, doublereal * alpha, doublereal * beta, doublereal * a, integer * lda, ftnlen uplo_len)
{
	/* System generated locals */
	integer a_dim1, a_offset, i__1, i__2, i__3;

	/* Local variables */
	static integer i__, j;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     October 31, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLASET initializes an m-by-n matrix A to BETA on the diagonal and */
	/*  ALPHA on the offdiagonals. */

	/*  Arguments */
	/*  ========= */

	/*  UPLO    (input) CHARACTER*1 */
	/*          Specifies the part of the matrix A to be set. */
	/*          = 'U':      Upper triangular part is set; the strictly lower */
	/*                      triangular part of A is not changed. */
	/*          = 'L':      Lower triangular part is set; the strictly upper */
	/*                      triangular part of A is not changed. */
	/*          Otherwise:  All of the matrix A is set. */

	/*  M       (input) INTEGER */
	/*          The number of rows of the matrix A.  M >= 0. */

	/*  N       (input) INTEGER */
	/*          The number of columns of the matrix A.  N >= 0. */

	/*  ALPHA   (input) DOUBLE PRECISION */
	/*          The constant to which the offdiagonal elements are to be set. */

	/*  BETA    (input) DOUBLE PRECISION */
	/*          The constant to which the diagonal elements are to be set. */

	/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          On exit, the leading m-by-n submatrix of A is set as follows: */

	/*          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n, */
	/*          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n, */
	/*          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j, */

	/*          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n). */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max(1,M). */

	/* ===================================================================== */

	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;

	/* Function Body */
	if (lsame_ (uplo, "U", (ftnlen) 1, (ftnlen) 1))
	{

		/*        Set the strictly upper triangular or trapezoidal part of the */
		/*        array to ALPHA. */

		i__1 = *n;
		for (j = 2; j <= i__1; ++j)
		{
			/* Computing MIN */
			i__3 = j - 1;
			i__2 = min (i__3, *m);
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				a[i__ + j * a_dim1] = *alpha;
				/* L10: */
			}
			/* L20: */
		}

	}
	else if (lsame_ (uplo, "L", (ftnlen) 1, (ftnlen) 1))
	{

		/*        Set the strictly lower triangular or trapezoidal part of the */
		/*        array to ALPHA. */

		i__1 = min (*m, *n);
		for (j = 1; j <= i__1; ++j)
		{
			i__2 = *m;
			for (i__ = j + 1; i__ <= i__2; ++i__)
			{
				a[i__ + j * a_dim1] = *alpha;
				/* L30: */
			}
			/* L40: */
		}

	}
	else
	{

		/*        Set the leading m-by-n submatrix to ALPHA. */

		i__1 = *n;
		for (j = 1; j <= i__1; ++j)
		{
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				a[i__ + j * a_dim1] = *alpha;
				/* L50: */
			}
			/* L60: */
		}
	}

	/*     Set the first min(M,N) diagonal elements to BETA. */

	i__1 = min (*m, *n);
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		a[i__ + i__ * a_dim1] = *beta;
		/* L70: */
	}

	return 0;

	/*     End of DLASET */

}	/* dlaset_ */

static int
dlassq_ (integer * n, doublereal * x, integer * incx, doublereal * scale, doublereal * sumsq)
{
	/* System generated locals */
	integer i__1, i__2;
	doublereal d__1;

	/* Local variables */
	static integer ix;
	static doublereal absxi;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     June 30, 1999 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLASSQ  returns the values  scl  and  smsq  such that */

	/*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq, */

	/*  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is */
	/*  assumed to be non-negative and  scl  returns the value */

	/*     scl = max( scale, abs( x( i ) ) ). */

	/*  scale and sumsq must be supplied in SCALE and SUMSQ and */
	/*  scl and smsq are overwritten on SCALE and SUMSQ respectively. */

	/*  The routine makes only one pass through the vector x. */

	/*  Arguments */
	/*  ========= */

	/*  N       (input) INTEGER */
	/*          The number of elements to be used from the vector X. */

	/*  X       (input) DOUBLE PRECISION array, dimension (N) */
	/*          The vector for which a scaled sum of squares is computed. */
	/*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n. */

	/*  INCX    (input) INTEGER */
	/*          The increment between successive values of the vector X. */
	/*          INCX > 0. */

	/*  SCALE   (input/output) DOUBLE PRECISION */
	/*          On entry, the value  scale  in the equation above. */
	/*          On exit, SCALE is overwritten with  scl , the scaling factor */
	/*          for the sum of squares. */

	/*  SUMSQ   (input/output) DOUBLE PRECISION */
	/*          On entry, the value  sumsq  in the equation above. */
	/*          On exit, SUMSQ is overwritten with  smsq , the basic sum of */
	/*          squares from which  scl  has been factored out. */

	/* ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/* Parameter adjustments */
	--x;

	/* Function Body */
	if (*n > 0)
	{
		i__1 = (*n - 1) * *incx + 1;
		i__2 = *incx;
		for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2)
		{
			if (x[ix] != 0.)
			{
				absxi = (d__1 = x[ix], abs (d__1));
				if (*scale < absxi)
				{
					/* Computing 2nd power */
					d__1 = *scale / absxi;
					*sumsq = *sumsq * (d__1 * d__1) + 1;
					*scale = absxi;
				}
				else
				{
					/* Computing 2nd power */
					d__1 = absxi / *scale;
					*sumsq += d__1 * d__1;
				}
			}
			/* L10: */
		}
	}
	return 0;

	/*     End of DLASSQ */

}	/* dlassq_ */

static int
dlaswp_ (integer * n, doublereal * a, integer * lda, integer * k1, integer * k2, integer * ipiv, integer * incx)
{
	/* System generated locals */
	integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

	/* Local variables */
	static integer i__, j, k, i1, i2, n32, ip, ix, ix0, inc;
	static doublereal temp;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     June 30, 1999 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLASWP performs a series of row interchanges on the matrix A. */
	/*  One row interchange is initiated for each of rows K1 through K2 of A. */

	/*  Arguments */
	/*  ========= */

	/*  N       (input) INTEGER */
	/*          The number of columns of the matrix A. */

	/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          On entry, the matrix of column dimension N to which the row */
	/*          interchanges will be applied. */
	/*          On exit, the permuted matrix. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A. */

	/*  K1      (input) INTEGER */
	/*          The first element of IPIV for which a row interchange will */
	/*          be done. */

	/*  K2      (input) INTEGER */
	/*          The last element of IPIV for which a row interchange will */
	/*          be done. */

	/*  IPIV    (input) INTEGER array, dimension (M*abs(INCX)) */
	/*          The vector of pivot indices.  Only the elements in positions */
	/*          K1 through K2 of IPIV are accessed. */
	/*          IPIV(K) = L implies rows K and L are to be interchanged. */

	/*  INCX    (input) INTEGER */
	/*          The increment between successive values of IPIV.  If IPIV */
	/*          is negative, the pivots are applied in reverse order. */

	/*  Further Details */
	/*  =============== */

	/*  Modified by */
	/*   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA */

	/* ===================================================================== */

	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Interchange row I with row IPIV(I) for each of rows K1 through K2. */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--ipiv;

	/* Function Body */
	if (*incx > 0)
	{
		ix0 = *k1;
		i1 = *k1;
		i2 = *k2;
		inc = 1;
	}
	else if (*incx < 0)
	{
		ix0 = (1 - *k2) * *incx + 1;
		i1 = *k2;
		i2 = *k1;
		inc = -1;
	}
	else
	{
		return 0;
	}

	n32 = *n / 32 << 5;
	if (n32 != 0)
	{
		i__1 = n32;
		for (j = 1; j <= i__1; j += 32)
		{
			ix = ix0;
			i__2 = i2;
			i__3 = inc;
			for (i__ = i1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3)
			{
				ip = ipiv[ix];
				if (ip != i__)
				{
					i__4 = j + 31;
					for (k = j; k <= i__4; ++k)
					{
						temp = a[i__ + k * a_dim1];
						a[i__ + k * a_dim1] = a[ip + k * a_dim1];
						a[ip + k * a_dim1] = temp;
						/* L10: */
					}
				}
				ix += *incx;
				/* L20: */
			}
			/* L30: */
		}
	}
	if (n32 != *n)
	{
		++n32;
		ix = ix0;
		i__1 = i2;
		i__3 = inc;
		for (i__ = i1; i__3 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__3)
		{
			ip = ipiv[ix];
			if (ip != i__)
			{
				i__2 = *n;
				for (k = n32; k <= i__2; ++k)
				{
					temp = a[i__ + k * a_dim1];
					a[i__ + k * a_dim1] = a[ip + k * a_dim1];
					a[ip + k * a_dim1] = temp;
					/* L40: */
				}
			}
			ix += *incx;
			/* L50: */
		}
	}

	return 0;

	/*     End of DLASWP */

}	/* dlaswp_ */

static int
dlasy2_ (logical * ltranl, logical * ltranr, integer * isgn,
			integer * n1, integer * n2, doublereal * tl, integer * ldtl, doublereal *
			tr, integer * ldtr, doublereal * b, integer * ldb, doublereal * scale, doublereal * x, integer * ldx, doublereal * xnorm, integer * info)
{
	/* Initialized data */

	static integer locu12[4] = { 3, 4, 1, 2 };
	static integer locl21[4] = { 2, 1, 4, 3 };
	static integer locu22[4] = { 4, 3, 2, 1 };
	static logical xswpiv[4] = { FALSE_, FALSE_, TRUE_, TRUE_ };
	static logical bswpiv[4] = { FALSE_, TRUE_, FALSE_, TRUE_ };

	/* System generated locals */
	integer b_dim1, b_offset, tl_dim1, tl_offset, tr_dim1, tr_offset, x_dim1, x_offset;
	doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;

	/* Local variables */
	static integer i__, j, k;
	static doublereal x2[2], l21, u11, u12;
	static integer ip, jp;
	static doublereal u22, t16[16] /* was [4][4] */ , gam, bet, eps, sgn, tmp[4], tau1, btmp[4], smin;
	static integer ipiv;
	static doublereal temp;
	static integer jpiv[4];
	static doublereal xmax;
	static integer ipsv, jpsv;
	static logical bswap;
	static logical xswap;
	static doublereal smlnum;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     October 31, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLASY2 solves for the N1 by N2 matrix X, 1 <= N1,N2 <= 2, in */

	/*         op(TL)*X + ISGN*X*op(TR) = SCALE*B, */

	/*  where TL is N1 by N1, TR is N2 by N2, B is N1 by N2, and ISGN = 1 or */
	/*  -1.  op(T) = T or T', where T' denotes the transpose of T. */

	/*  Arguments */
	/*  ========= */

	/*  LTRANL  (input) LOGICAL */
	/*          On entry, LTRANL specifies the op(TL): */
	/*             = .FALSE., op(TL) = TL, */
	/*             = .TRUE., op(TL) = TL'. */

	/*  LTRANR  (input) LOGICAL */
	/*          On entry, LTRANR specifies the op(TR): */
	/*            = .FALSE., op(TR) = TR, */
	/*            = .TRUE., op(TR) = TR'. */

	/*  ISGN    (input) INTEGER */
	/*          On entry, ISGN specifies the sign of the equation */
	/*          as described before. ISGN may only be 1 or -1. */

	/*  N1      (input) INTEGER */
	/*          On entry, N1 specifies the order of matrix TL. */
	/*          N1 may only be 0, 1 or 2. */

	/*  N2      (input) INTEGER */
	/*          On entry, N2 specifies the order of matrix TR. */
	/*          N2 may only be 0, 1 or 2. */

	/*  TL      (input) DOUBLE PRECISION array, dimension (LDTL,2) */
	/*          On entry, TL contains an N1 by N1 matrix. */

	/*  LDTL    (input) INTEGER */
	/*          The leading dimension of the matrix TL. LDTL >= max(1,N1). */

	/*  TR      (input) DOUBLE PRECISION array, dimension (LDTR,2) */
	/*          On entry, TR contains an N2 by N2 matrix. */

	/*  LDTR    (input) INTEGER */
	/*          The leading dimension of the matrix TR. LDTR >= max(1,N2). */

	/*  B       (input) DOUBLE PRECISION array, dimension (LDB,2) */
	/*          On entry, the N1 by N2 matrix B contains the right-hand */
	/*          side of the equation. */

	/*  LDB     (input) INTEGER */
	/*          The leading dimension of the matrix B. LDB >= max(1,N1). */

	/*  SCALE   (output) DOUBLE PRECISION */
	/*          On exit, SCALE contains the scale factor. SCALE is chosen */
	/*          less than or equal to 1 to prevent the solution overflowing. */

	/*  X       (output) DOUBLE PRECISION array, dimension (LDX,2) */
	/*          On exit, X contains the N1 by N2 solution. */

	/*  LDX     (input) INTEGER */
	/*          The leading dimension of the matrix X. LDX >= max(1,N1). */

	/*  XNORM   (output) DOUBLE PRECISION */
	/*          On exit, XNORM is the infinity-norm of the solution. */

	/*  INFO    (output) INTEGER */
	/*          On exit, INFO is set to */
	/*             0: successful exit. */
	/*             1: TL and TR have too close eigenvalues, so TL or */
	/*                TR is perturbed to get a nonsingular equation. */
	/*          NOTE: In the interests of speed, this routine does not */
	/*                check the inputs for errors. */

	/* ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Local Arrays .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Data statements .. */
	/* Parameter adjustments */
	tl_dim1 = *ldtl;
	tl_offset = 1 + tl_dim1;
	tl -= tl_offset;
	tr_dim1 = *ldtr;
	tr_offset = 1 + tr_dim1;
	tr -= tr_offset;
	b_dim1 = *ldb;
	b_offset = 1 + b_dim1;
	b -= b_offset;
	x_dim1 = *ldx;
	x_offset = 1 + x_dim1;
	x -= x_offset;

	/* Function Body */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Do not check the input parameters for errors */

	*info = 0;

	/*     Quick return if possible */

	if (*n1 == 0 || *n2 == 0)
	{
		return 0;
	}

	/*     Set constants to control overflow */

	eps = dlamch_ ("P", (ftnlen) 1);
	smlnum = dlamch_ ("S", (ftnlen) 1) / eps;
	sgn = (doublereal) (*isgn);

	k = *n1 + *n1 + *n2 - 2;
	switch (k)
	{
	case 1:
		goto L10;
	case 2:
		goto L20;
	case 3:
		goto L30;
	case 4:
		goto L50;
	}

	/*     1 by 1: TL11*X + SGN*X*TR11 = B11 */

 L10:
	tau1 = tl[tl_dim1 + 1] + sgn * tr[tr_dim1 + 1];
	bet = abs (tau1);
	if (bet <= smlnum)
	{
		tau1 = smlnum;
		bet = smlnum;
		*info = 1;
	}

	*scale = 1.;
	gam = (d__1 = b[b_dim1 + 1], abs (d__1));
	if (smlnum * gam > bet)
	{
		*scale = 1. / gam;
	}

	x[x_dim1 + 1] = b[b_dim1 + 1] * *scale / tau1;
	*xnorm = (d__1 = x[x_dim1 + 1], abs (d__1));
	return 0;

	/*     1 by 2: */
	/*     TL11*[X11 X12] + ISGN*[X11 X12]*op[TR11 TR12]  = [B11 B12] */
	/*                                       [TR21 TR22] */

 L20:

	/* Computing MAX */
	/* Computing MAX */
	d__7 = (d__1 = tl[tl_dim1 + 1], abs (d__1)), d__8 = (d__2 = tr[tr_dim1 + 1], abs (d__2)), d__7 = max (d__7, d__8), d__8 = (d__3 = tr[(tr_dim1 <<
																																													  1) + 1], abs (d__3)),
		d__7 = max (d__7, d__8), d__8 = (d__4 = tr[tr_dim1 + 2], abs (d__4)), d__7 = max (d__7, d__8), d__8 = (d__5 = tr[(tr_dim1 << 1) + 2], abs (d__5));
	d__6 = eps * max (d__7, d__8);
	smin = max (d__6, smlnum);
	tmp[0] = tl[tl_dim1 + 1] + sgn * tr[tr_dim1 + 1];
	tmp[3] = tl[tl_dim1 + 1] + sgn * tr[(tr_dim1 << 1) + 2];
	if (*ltranr)
	{
		tmp[1] = sgn * tr[tr_dim1 + 2];
		tmp[2] = sgn * tr[(tr_dim1 << 1) + 1];
	}
	else
	{
		tmp[1] = sgn * tr[(tr_dim1 << 1) + 1];
		tmp[2] = sgn * tr[tr_dim1 + 2];
	}
	btmp[0] = b[b_dim1 + 1];
	btmp[1] = b[(b_dim1 << 1) + 1];
	goto L40;

	/*     2 by 1: */
	/*          op[TL11 TL12]*[X11] + ISGN* [X11]*TR11  = [B11] */
	/*            [TL21 TL22] [X21]         [X21]         [B21] */

 L30:
	/* Computing MAX */
	/* Computing MAX */
	d__7 = (d__1 = tr[tr_dim1 + 1], abs (d__1)), d__8 = (d__2 = tl[tl_dim1 + 1], abs (d__2)), d__7 = max (d__7, d__8), d__8 = (d__3 = tl[(tl_dim1 <<
																																													  1) + 1], abs (d__3)),
		d__7 = max (d__7, d__8), d__8 = (d__4 = tl[tl_dim1 + 2], abs (d__4)), d__7 = max (d__7, d__8), d__8 = (d__5 = tl[(tl_dim1 << 1) + 2], abs (d__5));
	d__6 = eps * max (d__7, d__8);
	smin = max (d__6, smlnum);
	tmp[0] = tl[tl_dim1 + 1] + sgn * tr[tr_dim1 + 1];
	tmp[3] = tl[(tl_dim1 << 1) + 2] + sgn * tr[tr_dim1 + 1];
	if (*ltranl)
	{
		tmp[1] = tl[(tl_dim1 << 1) + 1];
		tmp[2] = tl[tl_dim1 + 2];
	}
	else
	{
		tmp[1] = tl[tl_dim1 + 2];
		tmp[2] = tl[(tl_dim1 << 1) + 1];
	}
	btmp[0] = b[b_dim1 + 1];
	btmp[1] = b[b_dim1 + 2];
 L40:

	/*     Solve 2 by 2 system using complete pivoting. */
	/*     Set pivots less than SMIN to SMIN. */

	ipiv = idamax_ (&c__4, tmp, &c__1);
	u11 = tmp[ipiv - 1];
	if (abs (u11) <= smin)
	{
		*info = 1;
		u11 = smin;
	}
	u12 = tmp[locu12[ipiv - 1] - 1];
	l21 = tmp[locl21[ipiv - 1] - 1] / u11;
	u22 = tmp[locu22[ipiv - 1] - 1] - u12 * l21;
	xswap = xswpiv[ipiv - 1];
	bswap = bswpiv[ipiv - 1];
	if (abs (u22) <= smin)
	{
		*info = 1;
		u22 = smin;
	}
	if (bswap)
	{
		temp = btmp[1];
		btmp[1] = btmp[0] - l21 * temp;
		btmp[0] = temp;
	}
	else
	{
		btmp[1] -= l21 * btmp[0];
	}
	*scale = 1.;
	if (smlnum * 2. * abs (btmp[1]) > abs (u22) || smlnum * 2. * abs (btmp[0]) > abs (u11))
	{
		/* Computing MAX */
		d__1 = abs (btmp[0]), d__2 = abs (btmp[1]);
		*scale = .5 / max (d__1, d__2);
		btmp[0] *= *scale;
		btmp[1] *= *scale;
	}
	x2[1] = btmp[1] / u22;
	x2[0] = btmp[0] / u11 - u12 / u11 * x2[1];
	if (xswap)
	{
		temp = x2[1];
		x2[1] = x2[0];
		x2[0] = temp;
	}
	x[x_dim1 + 1] = x2[0];
	if (*n1 == 1)
	{
		x[(x_dim1 << 1) + 1] = x2[1];
		*xnorm = (d__1 = x[x_dim1 + 1], abs (d__1)) + (d__2 = x[(x_dim1 << 1) + 1], abs (d__2));
	}
	else
	{
		x[x_dim1 + 2] = x2[1];
		/* Computing MAX */
		d__3 = (d__1 = x[x_dim1 + 1], abs (d__1)), d__4 = (d__2 = x[x_dim1 + 2], abs (d__2));
		*xnorm = max (d__3, d__4);
	}
	return 0;

	/*     2 by 2: */
	/*     op[TL11 TL12]*[X11 X12] +ISGN* [X11 X12]*op[TR11 TR12] = [B11 B12] */
	/*       [TL21 TL22] [X21 X22]        [X21 X22]   [TR21 TR22]   [B21 B22] */

	/*     Solve equivalent 4 by 4 system using complete pivoting. */
	/*     Set pivots less than SMIN to SMIN. */

 L50:
	/* Computing MAX */
	d__5 = (d__1 = tr[tr_dim1 + 1], abs (d__1)), d__6 = (d__2 = tr[(tr_dim1 <<
																						 1) + 1], abs (d__2)), d__5 = max (d__5, d__6), d__6 = (d__3 =
																																								  tr[tr_dim1 + 2], abs (d__3)), d__5 =
		max (d__5, d__6), d__6 = (d__4 = tr[(tr_dim1 << 1) + 2], abs (d__4));
	smin = max (d__5, d__6);
	/* Computing MAX */
	d__5 = smin, d__6 = (d__1 = tl[tl_dim1 + 1], abs (d__1)), d__5 = max (d__5,
																								 d__6), d__6 = (d__2 = tl[(tl_dim1 << 1) + 1], abs (d__2)), d__5 =
		max (d__5, d__6), d__6 = (d__3 = tl[tl_dim1 + 2], abs (d__3)), d__5 = max (d__5, d__6), d__6 = (d__4 = tl[(tl_dim1 << 1) + 2], abs (d__4));
	smin = max (d__5, d__6);
	/* Computing MAX */
	d__1 = eps * smin;
	smin = max (d__1, smlnum);
	btmp[0] = 0.;
	dcopy_ (&c__16, btmp, &c__0, t16, &c__1);
	t16[0] = tl[tl_dim1 + 1] + sgn * tr[tr_dim1 + 1];
	t16[5] = tl[(tl_dim1 << 1) + 2] + sgn * tr[tr_dim1 + 1];
	t16[10] = tl[tl_dim1 + 1] + sgn * tr[(tr_dim1 << 1) + 2];
	t16[15] = tl[(tl_dim1 << 1) + 2] + sgn * tr[(tr_dim1 << 1) + 2];
	if (*ltranl)
	{
		t16[4] = tl[tl_dim1 + 2];
		t16[1] = tl[(tl_dim1 << 1) + 1];
		t16[14] = tl[tl_dim1 + 2];
		t16[11] = tl[(tl_dim1 << 1) + 1];
	}
	else
	{
		t16[4] = tl[(tl_dim1 << 1) + 1];
		t16[1] = tl[tl_dim1 + 2];
		t16[14] = tl[(tl_dim1 << 1) + 1];
		t16[11] = tl[tl_dim1 + 2];
	}
	if (*ltranr)
	{
		t16[8] = sgn * tr[(tr_dim1 << 1) + 1];
		t16[13] = sgn * tr[(tr_dim1 << 1) + 1];
		t16[2] = sgn * tr[tr_dim1 + 2];
		t16[7] = sgn * tr[tr_dim1 + 2];
	}
	else
	{
		t16[8] = sgn * tr[tr_dim1 + 2];
		t16[13] = sgn * tr[tr_dim1 + 2];
		t16[2] = sgn * tr[(tr_dim1 << 1) + 1];
		t16[7] = sgn * tr[(tr_dim1 << 1) + 1];
	}
	btmp[0] = b[b_dim1 + 1];
	btmp[1] = b[b_dim1 + 2];
	btmp[2] = b[(b_dim1 << 1) + 1];
	btmp[3] = b[(b_dim1 << 1) + 2];

	/*     Perform elimination */

	for (i__ = 1; i__ <= 3; ++i__)
	{
		xmax = 0.;
		for (ip = i__; ip <= 4; ++ip)
		{
			for (jp = i__; jp <= 4; ++jp)
			{
				if ((d__1 = t16[ip + (jp << 2) - 5], abs (d__1)) >= xmax)
				{
					xmax = (d__1 = t16[ip + (jp << 2) - 5], abs (d__1));
					ipsv = ip;
					jpsv = jp;
				}
				/* L60: */
			}
			/* L70: */
		}
		if (ipsv != i__)
		{
			dswap_ (&c__4, &t16[ipsv - 1], &c__4, &t16[i__ - 1], &c__4);
			temp = btmp[i__ - 1];
			btmp[i__ - 1] = btmp[ipsv - 1];
			btmp[ipsv - 1] = temp;
		}
		if (jpsv != i__)
		{
			dswap_ (&c__4, &t16[(jpsv << 2) - 4], &c__1, &t16[(i__ << 2) - 4], &c__1);
		}
		jpiv[i__ - 1] = jpsv;
		if ((d__1 = t16[i__ + (i__ << 2) - 5], abs (d__1)) < smin)
		{
			*info = 1;
			t16[i__ + (i__ << 2) - 5] = smin;
		}
		for (j = i__ + 1; j <= 4; ++j)
		{
			t16[j + (i__ << 2) - 5] /= t16[i__ + (i__ << 2) - 5];
			btmp[j - 1] -= t16[j + (i__ << 2) - 5] * btmp[i__ - 1];
			for (k = i__ + 1; k <= 4; ++k)
			{
				t16[j + (k << 2) - 5] -= t16[j + (i__ << 2) - 5] * t16[i__ + (k << 2) - 5];
				/* L80: */
			}
			/* L90: */
		}
		/* L100: */
	}
	if (abs (t16[15]) < smin)
	{
		t16[15] = smin;
	}
	*scale = 1.;
	if (smlnum * 8. * abs (btmp[0]) > abs (t16[0]) || smlnum * 8. * abs (btmp[1])
		 > abs (t16[5]) || smlnum * 8. * abs (btmp[2]) > abs (t16[10]) || smlnum * 8. * abs (btmp[3]) > abs (t16[15]))
	{
		/* Computing MAX */
		d__1 = abs (btmp[0]), d__2 = abs (btmp[1]), d__1 = max (d__1, d__2), d__2 = abs (btmp[2]), d__1 = max (d__1, d__2), d__2 = abs (btmp[3]);
		*scale = .125 / max (d__1, d__2);
		btmp[0] *= *scale;
		btmp[1] *= *scale;
		btmp[2] *= *scale;
		btmp[3] *= *scale;
	}
	for (i__ = 1; i__ <= 4; ++i__)
	{
		k = 5 - i__;
		temp = 1. / t16[k + (k << 2) - 5];
		tmp[k - 1] = btmp[k - 1] * temp;
		for (j = k + 1; j <= 4; ++j)
		{
			tmp[k - 1] -= temp * t16[k + (j << 2) - 5] * tmp[j - 1];
			/* L110: */
		}
		/* L120: */
	}
	for (i__ = 1; i__ <= 3; ++i__)
	{
		if (jpiv[4 - i__ - 1] != 4 - i__)
		{
			temp = tmp[4 - i__ - 1];
			tmp[4 - i__ - 1] = tmp[jpiv[4 - i__ - 1] - 1];
			tmp[jpiv[4 - i__ - 1] - 1] = temp;
		}
		/* L130: */
	}
	x[x_dim1 + 1] = tmp[0];
	x[x_dim1 + 2] = tmp[1];
	x[(x_dim1 << 1) + 1] = tmp[2];
	x[(x_dim1 << 1) + 2] = tmp[3];
	/* Computing MAX */
	d__1 = abs (tmp[0]) + abs (tmp[2]), d__2 = abs (tmp[1]) + abs (tmp[3]);
	*xnorm = max (d__1, d__2);
	return 0;

	/*     End of DLASY2 */

}	/* dlasy2_ */

static int
dlatrs_ (char *uplo, char *trans, char *diag, char *normin, integer * n, doublereal * a, integer * lda, doublereal * x,
			doublereal * scale, doublereal * cnorm, integer * info, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len, ftnlen normin_len)
{
	/* System generated locals */
	integer a_dim1, a_offset, i__1, i__2, i__3;
	doublereal d__1, d__2, d__3;

	/* Local variables */
	static integer i__, j;
	static doublereal xj, rec, tjj;
	static integer jinc;
	static doublereal xbnd;
	static integer imax;
	static doublereal tmax, tjjs, xmax, grow, sumj;
	static doublereal tscal, uscal;
	static integer jlast;
	static logical upper;
	static doublereal bignum;
	static logical notran;
	static integer jfirst;
	static doublereal smlnum;
	static logical nounit;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     June 30, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DLATRS solves one of the triangular systems */

	/*     A *x = s*b  or  A'*x = s*b */

	/*  with scaling to prevent overflow.  Here A is an upper or lower */
	/*  triangular matrix, A' denotes the transpose of A, x and b are */
	/*  n-element vectors, and s is a scaling factor, usually less than */
	/*  or equal to 1, chosen so that the components of x will be less than */
	/*  the overflow threshold.  If the unscaled problem will not cause */
	/*  overflow, the Level 2 BLAS routine DTRSV is called.  If the matrix A */
	/*  is singular (A(j,j) = 0 for some j), then s is set to 0 and a */
	/*  non-trivial solution to A*x = 0 is returned. */

	/*  Arguments */
	/*  ========= */

	/*  UPLO    (input) CHARACTER*1 */
	/*          Specifies whether the matrix A is upper or lower triangular. */
	/*          = 'U':  Upper triangular */
	/*          = 'L':  Lower triangular */

	/*  TRANS   (input) CHARACTER*1 */
	/*          Specifies the operation applied to A. */
	/*          = 'N':  Solve A * x = s*b  (No transpose) */
	/*          = 'T':  Solve A'* x = s*b  (Transpose) */
	/*          = 'C':  Solve A'* x = s*b  (Conjugate transpose = Transpose) */

	/*  DIAG    (input) CHARACTER*1 */
	/*          Specifies whether or not the matrix A is unit triangular. */
	/*          = 'N':  Non-unit triangular */
	/*          = 'U':  Unit triangular */

	/*  NORMIN  (input) CHARACTER*1 */
	/*          Specifies whether CNORM has been set or not. */
	/*          = 'Y':  CNORM contains the column norms on entry */
	/*          = 'N':  CNORM is not set on entry.  On exit, the norms will */
	/*                  be computed and stored in CNORM. */

	/*  N       (input) INTEGER */
	/*          The order of the matrix A.  N >= 0. */

	/*  A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          The triangular matrix A.  If UPLO = 'U', the leading n by n */
	/*          upper triangular part of the array A contains the upper */
	/*          triangular matrix, and the strictly lower triangular part of */
	/*          A is not referenced.  If UPLO = 'L', the leading n by n lower */
	/*          triangular part of the array A contains the lower triangular */
	/*          matrix, and the strictly upper triangular part of A is not */
	/*          referenced.  If DIAG = 'U', the diagonal elements of A are */
	/*          also not referenced and are assumed to be 1. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A.  LDA >= max (1,N). */

	/*  X       (input/output) DOUBLE PRECISION array, dimension (N) */
	/*          On entry, the right hand side b of the triangular system. */
	/*          On exit, X is overwritten by the solution vector x. */

	/*  SCALE   (output) DOUBLE PRECISION */
	/*          The scaling factor s for the triangular system */
	/*             A * x = s*b  or  A'* x = s*b. */
	/*          If SCALE = 0, the matrix A is singular or badly scaled, and */
	/*          the vector x is an exact or approximate solution to A*x = 0. */

	/*  CNORM   (input or output) DOUBLE PRECISION array, dimension (N) */

	/*          If NORMIN = 'Y', CNORM is an input argument and CNORM(j) */
	/*          contains the norm of the off-diagonal part of the j-th column */
	/*          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal */
	/*          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j) */
	/*          must be greater than or equal to the 1-norm. */

	/*          If NORMIN = 'N', CNORM is an output argument and CNORM(j) */
	/*          returns the 1-norm of the offdiagonal part of the j-th column */
	/*          of A. */

	/*  INFO    (output) INTEGER */
	/*          = 0:  successful exit */
	/*          < 0:  if INFO = -k, the k-th argument had an illegal value */

	/*  Further Details */
	/*  ======= ======= */

	/*  A rough bound on x is computed; if that is less than overflow, DTRSV */
	/*  is called, otherwise, specific code is used which checks for possible */
	/*  overflow or divide-by-zero at every operation. */

	/*  A columnwise scheme is used for solving A*x = b.  The basic algorithm */
	/*  if A is lower triangular is */

	/*       x[1:n] := b[1:n] */
	/*       for j = 1, ..., n */
	/*            x(j) := x(j) / A(j,j) */
	/*            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j] */
	/*       end */

	/*  Define bounds on the components of x after j iterations of the loop: */
	/*     M(j) = bound on x[1:j] */
	/*     G(j) = bound on x[j+1:n] */
	/*  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}. */

	/*  Then for iteration j+1 we have */
	/*     M(j+1) <= G(j) / | A(j+1,j+1) | */
	/*     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] | */
	/*            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | ) */

	/*  where CNORM(j+1) is greater than or equal to the infinity-norm of */
	/*  column j+1 of A, not counting the diagonal.  Hence */

	/*     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | ) */
	/*                  1<=i<=j */
	/*  and */

	/*     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| ) */
	/*                                   1<=i< j */

	/*  Since |x(j)| <= M(j), we use the Level 2 BLAS routine DTRSV if the */
	/*  reciprocal of the largest M(j), j=1,..,n, is larger than */
	/*  max(underflow, 1/overflow). */

	/*  The bound on x(j) is also used to determine when a step in the */
	/*  columnwise method can be performed without fear of overflow.  If */
	/*  the computed bound is greater than a large constant, x is scaled to */
	/*  prevent overflow, but if the bound overflows, x is set to 0, x(j) to */
	/*  1, and scale to 0, and a non-trivial solution to A*x = 0 is found. */

	/*  Similarly, a row-wise scheme is used to solve A'*x = b.  The basic */
	/*  algorithm for A upper triangular is */

	/*       for j = 1, ..., n */
	/*            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j) */
	/*       end */

	/*  We simultaneously compute two bounds */
	/*       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j */
	/*       M(j) = bound on x(i), 1<=i<=j */

	/*  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we */
	/*  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1. */
	/*  Then the bound on x(j) is */

	/*       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) | */

	/*            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| ) */
	/*                      1<=i<=j */

	/*  and we can safely call DTRSV if 1/M(n) and 1/G(n) are both greater */
	/*  than max(underflow, 1/overflow). */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--x;
	--cnorm;

	/* Function Body */
	*info = 0;
	upper = lsame_ (uplo, "U", (ftnlen) 1, (ftnlen) 1);
	notran = lsame_ (trans, "N", (ftnlen) 1, (ftnlen) 1);
	nounit = lsame_ (diag, "N", (ftnlen) 1, (ftnlen) 1);

	/*     Test the input parameters. */

	if (!upper && !lsame_ (uplo, "L", (ftnlen) 1, (ftnlen) 1))
	{
		*info = -1;
	}
	else if (!notran && !lsame_ (trans, "T", (ftnlen) 1, (ftnlen) 1) && !lsame_ (trans, "C", (ftnlen) 1, (ftnlen) 1))
	{
		*info = -2;
	}
	else if (!nounit && !lsame_ (diag, "U", (ftnlen) 1, (ftnlen) 1))
	{
		*info = -3;
	}
	else if (!lsame_ (normin, "Y", (ftnlen) 1, (ftnlen) 1) && !lsame_ (normin, "N", (ftnlen) 1, (ftnlen) 1))
	{
		*info = -4;
	}
	else if (*n < 0)
	{
		*info = -5;
	}
	else if (*lda < max (1, *n))
	{
		*info = -7;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DLATRS", &i__1, (ftnlen) 6);
		return 0;
	}

	/*     Quick return if possible */

	if (*n == 0)
	{
		return 0;
	}

	/*     Determine machine dependent parameters to control overflow. */

	smlnum = dlamch_ ("Safe minimum", (ftnlen) 12) / dlamch_ ("Precision", (ftnlen) 9);
	bignum = 1. / smlnum;
	*scale = 1.;

	if (lsame_ (normin, "N", (ftnlen) 1, (ftnlen) 1))
	{

		/*        Compute the 1-norm of each column, not including the diagonal. */

		if (upper)
		{

			/*           A is upper triangular. */

			i__1 = *n;
			for (j = 1; j <= i__1; ++j)
			{
				i__2 = j - 1;
				cnorm[j] = dasum_ (&i__2, &a[j * a_dim1 + 1], &c__1);
				/* L10: */
			}
		}
		else
		{

			/*           A is lower triangular. */

			i__1 = *n - 1;
			for (j = 1; j <= i__1; ++j)
			{
				i__2 = *n - j;
				cnorm[j] = dasum_ (&i__2, &a[j + 1 + j * a_dim1], &c__1);
				/* L20: */
			}
			cnorm[*n] = 0.;
		}
	}

	/*     Scale the column norms by TSCAL if the maximum element in CNORM is */
	/*     greater than BIGNUM. */

	imax = idamax_ (n, &cnorm[1], &c__1);
	tmax = cnorm[imax];
	if (tmax <= bignum)
	{
		tscal = 1.;
	}
	else
	{
		tscal = 1. / (smlnum * tmax);
		dscal_ (n, &tscal, &cnorm[1], &c__1);
	}

	/*     Compute a bound on the computed solution vector to see if the */
	/*     Level 2 BLAS routine DTRSV can be used. */

	j = idamax_ (n, &x[1], &c__1);
	xmax = (d__1 = x[j], abs (d__1));
	xbnd = xmax;
	if (notran)
	{

		/*        Compute the growth in A * x = b. */

		if (upper)
		{
			jfirst = *n;
			jlast = 1;
			jinc = -1;
		}
		else
		{
			jfirst = 1;
			jlast = *n;
			jinc = 1;
		}

		if (tscal != 1.)
		{
			grow = 0.;
			goto L50;
		}

		if (nounit)
		{

			/*           A is non-unit triangular. */

			/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
			/*           Initially, G(0) = max{x(i), i=1,...,n}. */

			grow = 1. / max (xbnd, smlnum);
			xbnd = grow;
			i__1 = jlast;
			i__2 = jinc;
			for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2)
			{

				/*              Exit the loop if the growth factor is too small. */

				if (grow <= smlnum)
				{
					goto L50;
				}

				/*              M(j) = G(j-1) / abs(A(j,j)) */

				tjj = (d__1 = a[j + j * a_dim1], abs (d__1));
				/* Computing MIN */
				d__1 = xbnd, d__2 = min (1., tjj) * grow;
				xbnd = min (d__1, d__2);
				if (tjj + cnorm[j] >= smlnum)
				{

					/*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) ) */

					grow *= tjj / (tjj + cnorm[j]);
				}
				else
				{

					/*                 G(j) could overflow, set GROW to 0. */

					grow = 0.;
				}
				/* L30: */
			}
			grow = xbnd;
		}
		else
		{

			/*           A is unit triangular. */

			/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

			/* Computing MIN */
			d__1 = 1., d__2 = 1. / max (xbnd, smlnum);
			grow = min (d__1, d__2);
			i__2 = jlast;
			i__1 = jinc;
			for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1)
			{

				/*              Exit the loop if the growth factor is too small. */

				if (grow <= smlnum)
				{
					goto L50;
				}

				/*              G(j) = G(j-1)*( 1 + CNORM(j) ) */

				grow *= 1. / (cnorm[j] + 1.);
				/* L40: */
			}
		}
	 L50:

		;
	}
	else
	{

		/*        Compute the growth in A' * x = b. */

		if (upper)
		{
			jfirst = 1;
			jlast = *n;
			jinc = 1;
		}
		else
		{
			jfirst = *n;
			jlast = 1;
			jinc = -1;
		}

		if (tscal != 1.)
		{
			grow = 0.;
			goto L80;
		}

		if (nounit)
		{

			/*           A is non-unit triangular. */

			/*           Compute GROW = 1/G(j) and XBND = 1/M(j). */
			/*           Initially, M(0) = max{x(i), i=1,...,n}. */

			grow = 1. / max (xbnd, smlnum);
			xbnd = grow;
			i__1 = jlast;
			i__2 = jinc;
			for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2)
			{

				/*              Exit the loop if the growth factor is too small. */

				if (grow <= smlnum)
				{
					goto L80;
				}

				/*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) ) */

				xj = cnorm[j] + 1.;
				/* Computing MIN */
				d__1 = grow, d__2 = xbnd / xj;
				grow = min (d__1, d__2);

				/*              M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j)) */

				tjj = (d__1 = a[j + j * a_dim1], abs (d__1));
				if (xj > tjj)
				{
					xbnd *= tjj / xj;
				}
				/* L60: */
			}
			grow = min (grow, xbnd);
		}
		else
		{

			/*           A is unit triangular. */

			/*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}. */

			/* Computing MIN */
			d__1 = 1., d__2 = 1. / max (xbnd, smlnum);
			grow = min (d__1, d__2);
			i__2 = jlast;
			i__1 = jinc;
			for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1)
			{

				/*              Exit the loop if the growth factor is too small. */

				if (grow <= smlnum)
				{
					goto L80;
				}

				/*              G(j) = ( 1 + CNORM(j) )*G(j-1) */

				xj = cnorm[j] + 1.;
				grow /= xj;
				/* L70: */
			}
		}
	 L80:
		;
	}

	if (grow * tscal > smlnum)
	{

		/*        Use the Level 2 BLAS solve if the reciprocal of the bound on */
		/*        elements of X is not too small. */

		dtrsv_ (uplo, trans, diag, n, &a[a_offset], lda, &x[1], &c__1, (ftnlen) 1, (ftnlen) 1, (ftnlen) 1);
	}
	else
	{

		/*        Use a Level 1 BLAS solve, scaling intermediate results. */

		if (xmax > bignum)
		{

			/*           Scale X so that its components are less than or equal to */
			/*           BIGNUM in absolute value. */

			*scale = bignum / xmax;
			dscal_ (n, scale, &x[1], &c__1);
			xmax = bignum;
		}

		if (notran)
		{

			/*           Solve A * x = b */

			i__1 = jlast;
			i__2 = jinc;
			for (j = jfirst; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2)
			{

				/*              Compute x(j) = b(j) / A(j,j), scaling x if necessary. */

				xj = (d__1 = x[j], abs (d__1));
				if (nounit)
				{
					tjjs = a[j + j * a_dim1] * tscal;
				}
				else
				{
					tjjs = tscal;
					if (tscal == 1.)
					{
						goto L100;
					}
				}
				tjj = abs (tjjs);
				if (tjj > smlnum)
				{

					/*                    abs(A(j,j)) > SMLNUM: */

					if (tjj < 1.)
					{
						if (xj > tjj * bignum)
						{

							/*                          Scale x by 1/b(j). */

							rec = 1. / xj;
							dscal_ (n, &rec, &x[1], &c__1);
							*scale *= rec;
							xmax *= rec;
						}
					}
					x[j] /= tjjs;
					xj = (d__1 = x[j], abs (d__1));
				}
				else if (tjj > 0.)
				{

					/*                    0 < abs(A(j,j)) <= SMLNUM: */

					if (xj > tjj * bignum)
					{

						/*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM */
						/*                       to avoid overflow when dividing by A(j,j). */

						rec = tjj * bignum / xj;
						if (cnorm[j] > 1.)
						{

							/*                          Scale by 1/CNORM(j) to avoid overflow when */
							/*                          multiplying x(j) times column j. */

							rec /= cnorm[j];
						}
						dscal_ (n, &rec, &x[1], &c__1);
						*scale *= rec;
						xmax *= rec;
					}
					x[j] /= tjjs;
					xj = (d__1 = x[j], abs (d__1));
				}
				else
				{

					/*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
					/*                    scale = 0, and compute a solution to A*x = 0. */

					i__3 = *n;
					for (i__ = 1; i__ <= i__3; ++i__)
					{
						x[i__] = 0.;
						/* L90: */
					}
					x[j] = 1.;
					xj = 1.;
					*scale = 0.;
					xmax = 0.;
				}
			 L100:

				/*              Scale x if necessary to avoid overflow when adding a */
				/*              multiple of column j of A. */

				if (xj > 1.)
				{
					rec = 1. / xj;
					if (cnorm[j] > (bignum - xmax) * rec)
					{

						/*                    Scale x by 1/(2*abs(x(j))). */

						rec *= .5;
						dscal_ (n, &rec, &x[1], &c__1);
						*scale *= rec;
					}
				}
				else if (xj * cnorm[j] > bignum - xmax)
				{

					/*                 Scale x by 1/2. */

					dscal_ (n, &c_b1496, &x[1], &c__1);
					*scale *= .5;
				}

				if (upper)
				{
					if (j > 1)
					{

						/*                    Compute the update */
						/*                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j) */

						i__3 = j - 1;
						d__1 = -x[j] * tscal;
						daxpy_ (&i__3, &d__1, &a[j * a_dim1 + 1], &c__1, &x[1], &c__1);
						i__3 = j - 1;
						i__ = idamax_ (&i__3, &x[1], &c__1);
						xmax = (d__1 = x[i__], abs (d__1));
					}
				}
				else
				{
					if (j < *n)
					{

						/*                    Compute the update */
						/*                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j) */

						i__3 = *n - j;
						d__1 = -x[j] * tscal;
						daxpy_ (&i__3, &d__1, &a[j + 1 + j * a_dim1], &c__1, &x[j + 1], &c__1);
						i__3 = *n - j;
						i__ = j + idamax_ (&i__3, &x[j + 1], &c__1);
						xmax = (d__1 = x[i__], abs (d__1));
					}
				}
				/* L110: */
			}

		}
		else
		{

			/*           Solve A' * x = b */

			i__2 = jlast;
			i__1 = jinc;
			for (j = jfirst; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1)
			{

				/*              Compute x(j) = b(j) - sum A(k,j)*x(k). */
				/*                                    k<>j */

				xj = (d__1 = x[j], abs (d__1));
				uscal = tscal;
				rec = 1. / max (xmax, 1.);
				if (cnorm[j] > (bignum - xj) * rec)
				{

					/*                 If x(j) could overflow, scale x by 1/(2*XMAX). */

					rec *= .5;
					if (nounit)
					{
						tjjs = a[j + j * a_dim1] * tscal;
					}
					else
					{
						tjjs = tscal;
					}
					tjj = abs (tjjs);
					if (tjj > 1.)
					{

						/*                       Divide by A(j,j) when scaling x if A(j,j) > 1. */

						/* Computing MIN */
						d__1 = 1., d__2 = rec * tjj;
						rec = min (d__1, d__2);
						uscal /= tjjs;
					}
					if (rec < 1.)
					{
						dscal_ (n, &rec, &x[1], &c__1);
						*scale *= rec;
						xmax *= rec;
					}
				}

				sumj = 0.;
				if (uscal == 1.)
				{

					/*                 If the scaling needed for A in the dot product is 1, */
					/*                 call DDOT to perform the dot product. */

					if (upper)
					{
						i__3 = j - 1;
						sumj = ddot_ (&i__3, &a[j * a_dim1 + 1], &c__1, &x[1], &c__1);
					}
					else if (j < *n)
					{
						i__3 = *n - j;
						sumj = ddot_ (&i__3, &a[j + 1 + j * a_dim1], &c__1, &x[j + 1], &c__1);
					}
				}
				else
				{

					/*                 Otherwise, use in-line code for the dot product. */

					if (upper)
					{
						i__3 = j - 1;
						for (i__ = 1; i__ <= i__3; ++i__)
						{
							sumj += a[i__ + j * a_dim1] * uscal * x[i__];
							/* L120: */
						}
					}
					else if (j < *n)
					{
						i__3 = *n;
						for (i__ = j + 1; i__ <= i__3; ++i__)
						{
							sumj += a[i__ + j * a_dim1] * uscal * x[i__];
							/* L130: */
						}
					}
				}

				if (uscal == tscal)
				{

					/*                 Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j) */
					/*                 was not used to scale the dotproduct. */

					x[j] -= sumj;
					xj = (d__1 = x[j], abs (d__1));
					if (nounit)
					{
						tjjs = a[j + j * a_dim1] * tscal;
					}
					else
					{
						tjjs = tscal;
						if (tscal == 1.)
						{
							goto L150;
						}
					}

					/*                    Compute x(j) = x(j) / A(j,j), scaling if necessary. */

					tjj = abs (tjjs);
					if (tjj > smlnum)
					{

						/*                       abs(A(j,j)) > SMLNUM: */

						if (tjj < 1.)
						{
							if (xj > tjj * bignum)
							{

								/*                             Scale X by 1/abs(x(j)). */

								rec = 1. / xj;
								dscal_ (n, &rec, &x[1], &c__1);
								*scale *= rec;
								xmax *= rec;
							}
						}
						x[j] /= tjjs;
					}
					else if (tjj > 0.)
					{

						/*                       0 < abs(A(j,j)) <= SMLNUM: */

						if (xj > tjj * bignum)
						{

							/*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM. */

							rec = tjj * bignum / xj;
							dscal_ (n, &rec, &x[1], &c__1);
							*scale *= rec;
							xmax *= rec;
						}
						x[j] /= tjjs;
					}
					else
					{

						/*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and */
						/*                       scale = 0, and compute a solution to A'*x = 0. */

						i__3 = *n;
						for (i__ = 1; i__ <= i__3; ++i__)
						{
							x[i__] = 0.;
							/* L140: */
						}
						x[j] = 1.;
						*scale = 0.;
						xmax = 0.;
					}
				 L150:
					;
				}
				else
				{

					/*                 Compute x(j) := x(j) / A(j,j)  - sumj if the dot */
					/*                 product has already been divided by 1/A(j,j). */

					x[j] = x[j] / tjjs - sumj;
				}
				/* Computing MAX */
				d__2 = xmax, d__3 = (d__1 = x[j], abs (d__1));
				xmax = max (d__2, d__3);
				/* L160: */
			}
		}
		*scale /= tscal;
	}

	/*     Scale the column norms by 1/TSCAL for return. */

	if (tscal != 1.)
	{
		d__1 = 1. / tscal;
		dscal_ (n, &d__1, &cnorm[1], &c__1);
	}

	return 0;

	/*     End of DLATRS */

}	/* dlatrs_ */

static int
dorg2r_ (integer * m, integer * n, integer * k, doublereal * a, integer * lda, doublereal * tau, doublereal * work, integer * info)
{
	/* System generated locals */
	integer a_dim1, a_offset, i__1, i__2;
	doublereal d__1;

	/* Local variables */
	static integer i__, j, l;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     February 29, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DORG2R generates an m by n real matrix Q with orthonormal columns, */
	/*  which is defined as the first n columns of a product of k elementary */
	/*  reflectors of order m */

	/*        Q  =  H(1) H(2) . . . H(k) */

	/*  as returned by DGEQRF. */

	/*  Arguments */
	/*  ========= */

	/*  M       (input) INTEGER */
	/*          The number of rows of the matrix Q. M >= 0. */

	/*  N       (input) INTEGER */
	/*          The number of columns of the matrix Q. M >= N >= 0. */

	/*  K       (input) INTEGER */
	/*          The number of elementary reflectors whose product defines the */
	/*          matrix Q. N >= K >= 0. */

	/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          On entry, the i-th column must contain the vector which */
	/*          defines the elementary reflector H(i), for i = 1,2,...,k, as */
	/*          returned by DGEQRF in the first k columns of its array */
	/*          argument A. */
	/*          On exit, the m-by-n matrix Q. */

	/*  LDA     (input) INTEGER */
	/*          The first dimension of the array A. LDA >= max(1,M). */

	/*  TAU     (input) DOUBLE PRECISION array, dimension (K) */
	/*          TAU(i) must contain the scalar factor of the elementary */
	/*          reflector H(i), as returned by DGEQRF. */

	/*  WORK    (workspace) DOUBLE PRECISION array, dimension (N) */

	/*  INFO    (output) INTEGER */
	/*          = 0: successful exit */
	/*          < 0: if INFO = -i, the i-th argument has an illegal value */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Test the input arguments */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--tau;
	--work;

	/* Function Body */
	*info = 0;
	if (*m < 0)
	{
		*info = -1;
	}
	else if (*n < 0 || *n > *m)
	{
		*info = -2;
	}
	else if (*k < 0 || *k > *n)
	{
		*info = -3;
	}
	else if (*lda < max (1, *m))
	{
		*info = -5;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DORG2R", &i__1, (ftnlen) 6);
		return 0;
	}

	/*     Quick return if possible */

	if (*n <= 0)
	{
		return 0;
	}

	/*     Initialise columns k+1:n to columns of the unit matrix */

	i__1 = *n;
	for (j = *k + 1; j <= i__1; ++j)
	{
		i__2 = *m;
		for (l = 1; l <= i__2; ++l)
		{
			a[l + j * a_dim1] = 0.;
			/* L10: */
		}
		a[j + j * a_dim1] = 1.;
		/* L20: */
	}

	for (i__ = *k; i__ >= 1; --i__)
	{

		/*        Apply H(i) to A(i:m,i:n) from the left */

		if (i__ < *n)
		{
			a[i__ + i__ * a_dim1] = 1.;
			i__1 = *m - i__ + 1;
			i__2 = *n - i__;
			dlarf_ ("Left", &i__1, &i__2, &a[i__ + i__ * a_dim1], &c__1, &tau[i__], &a[i__ + (i__ + 1) * a_dim1], lda, &work[1], (ftnlen) 4);
		}
		if (i__ < *m)
		{
			i__1 = *m - i__;
			d__1 = -tau[i__];
			dscal_ (&i__1, &d__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
		}
		a[i__ + i__ * a_dim1] = 1. - tau[i__];

		/*        Set A(1:i-1,i) to zero */

		i__1 = i__ - 1;
		for (l = 1; l <= i__1; ++l)
		{
			a[l + i__ * a_dim1] = 0.;
			/* L30: */
		}
		/* L40: */
	}
	return 0;

	/*     End of DORG2R */

}	/* dorg2r_ */

static int
dorghr_ (integer * n, integer * ilo, integer * ihi, doublereal * a, integer * lda, doublereal * tau, doublereal * work, integer * lwork, integer * info)
{
	/* System generated locals */
	integer a_dim1, a_offset, i__1, i__2;

	/* Local variables */
	static integer i__, j, nb, nh, iinfo;
	static integer lwkopt;
	static logical lquery;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     June 30, 1999 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DORGHR generates a real orthogonal matrix Q which is defined as the */
	/*  product of IHI-ILO elementary reflectors of order N, as returned by */
	/*  DGEHRD: */

	/*  Q = H(ilo) H(ilo+1) . . . H(ihi-1). */

	/*  Arguments */
	/*  ========= */

	/*  N       (input) INTEGER */
	/*          The order of the matrix Q. N >= 0. */

	/*  ILO     (input) INTEGER */
	/*  IHI     (input) INTEGER */
	/*          ILO and IHI must have the same values as in the previous call */
	/*          of DGEHRD. Q is equal to the unit matrix except in the */
	/*          submatrix Q(ilo+1:ihi,ilo+1:ihi). */
	/*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0. */

	/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          On entry, the vectors which define the elementary reflectors, */
	/*          as returned by DGEHRD. */
	/*          On exit, the N-by-N orthogonal matrix Q. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A. LDA >= max(1,N). */

	/*  TAU     (input) DOUBLE PRECISION array, dimension (N-1) */
	/*          TAU(i) must contain the scalar factor of the elementary */
	/*          reflector H(i), as returned by DGEHRD. */

	/*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) */
	/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

	/*  LWORK   (input) INTEGER */
	/*          The dimension of the array WORK. LWORK >= IHI-ILO. */
	/*          For optimum performance LWORK >= (IHI-ILO)*NB, where NB is */
	/*          the optimal blocksize. */

	/*          If LWORK = -1, then a workspace query is assumed; the routine */
	/*          only calculates the optimal size of the WORK array, returns */
	/*          this value as the first entry of the WORK array, and no error */
	/*          message related to LWORK is issued by XERBLA. */

	/*  INFO    (output) INTEGER */
	/*          = 0:  successful exit */
	/*          < 0:  if INFO = -i, the i-th argument had an illegal value */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Test the input arguments */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--tau;
	--work;

	/* Function Body */
	*info = 0;
	nh = *ihi - *ilo;
	lquery = *lwork == -1;
	if (*n < 0)
	{
		*info = -1;
	}
	else if (*ilo < 1 || *ilo > max (1, *n))
	{
		*info = -2;
	}
	else if (*ihi < min (*ilo, *n) || *ihi > *n)
	{
		*info = -3;
	}
	else if (*lda < max (1, *n))
	{
		*info = -5;
	}
	else if (*lwork < max (1, nh) && !lquery)
	{
		*info = -8;
	}

	if (*info == 0)
	{
		nb = ilaenv_ (&c__1, "DORGQR", " ", &nh, &nh, &nh, &c_n1, (ftnlen) 6, (ftnlen) 1);
		lwkopt = max (1, nh) * nb;
		work[1] = (doublereal) lwkopt;
	}

	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DORGHR", &i__1, (ftnlen) 6);
		return 0;
	}
	else if (lquery)
	{
		return 0;
	}

	/*     Quick return if possible */

	if (*n == 0)
	{
		work[1] = 1.;
		return 0;
	}

	/*     Shift the vectors which define the elementary reflectors one */
	/*     column to the right, and set the first ilo and the last n-ihi */
	/*     rows and columns to those of the unit matrix */

	i__1 = *ilo + 1;
	for (j = *ihi; j >= i__1; --j)
	{
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__)
		{
			a[i__ + j * a_dim1] = 0.;
			/* L10: */
		}
		i__2 = *ihi;
		for (i__ = j + 1; i__ <= i__2; ++i__)
		{
			a[i__ + j * a_dim1] = a[i__ + (j - 1) * a_dim1];
			/* L20: */
		}
		i__2 = *n;
		for (i__ = *ihi + 1; i__ <= i__2; ++i__)
		{
			a[i__ + j * a_dim1] = 0.;
			/* L30: */
		}
		/* L40: */
	}
	i__1 = *ilo;
	for (j = 1; j <= i__1; ++j)
	{
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__)
		{
			a[i__ + j * a_dim1] = 0.;
			/* L50: */
		}
		a[j + j * a_dim1] = 1.;
		/* L60: */
	}
	i__1 = *n;
	for (j = *ihi + 1; j <= i__1; ++j)
	{
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__)
		{
			a[i__ + j * a_dim1] = 0.;
			/* L70: */
		}
		a[j + j * a_dim1] = 1.;
		/* L80: */
	}

	if (nh > 0)
	{

		/*        Generate Q(ilo+1:ihi,ilo+1:ihi) */

		dorgqr_ (&nh, &nh, &nh, &a[*ilo + 1 + (*ilo + 1) * a_dim1], lda, &tau[*ilo], &work[1], lwork, &iinfo);
	}
	work[1] = (doublereal) lwkopt;
	return 0;

	/*     End of DORGHR */

}	/* dorghr_ */

static int
dorgqr_ (integer * m, integer * n, integer * k, doublereal * a, integer * lda, doublereal * tau, doublereal * work, integer * lwork, integer * info)
{
	/* System generated locals */
	integer a_dim1, a_offset, i__1, i__2, i__3;

	/* Local variables */
	static integer i__, j, l, ib, nb, ki, kk, nx, iws, nbmin, iinfo;
	static integer ldwork, lwkopt;
	static logical lquery;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     June 30, 1999 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DORGQR generates an M-by-N real matrix Q with orthonormal columns, */
	/*  which is defined as the first N columns of a product of K elementary */
	/*  reflectors of order M */

	/*        Q  =  H(1) H(2) . . . H(k) */

	/*  as returned by DGEQRF. */

	/*  Arguments */
	/*  ========= */

	/*  M       (input) INTEGER */
	/*          The number of rows of the matrix Q. M >= 0. */

	/*  N       (input) INTEGER */
	/*          The number of columns of the matrix Q. M >= N >= 0. */

	/*  K       (input) INTEGER */
	/*          The number of elementary reflectors whose product defines the */
	/*          matrix Q. N >= K >= 0. */

	/*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
	/*          On entry, the i-th column must contain the vector which */
	/*          defines the elementary reflector H(i), for i = 1,2,...,k, as */
	/*          returned by DGEQRF in the first k columns of its array */
	/*          argument A. */
	/*          On exit, the M-by-N matrix Q. */

	/*  LDA     (input) INTEGER */
	/*          The first dimension of the array A. LDA >= max(1,M). */

	/*  TAU     (input) DOUBLE PRECISION array, dimension (K) */
	/*          TAU(i) must contain the scalar factor of the elementary */
	/*          reflector H(i), as returned by DGEQRF. */

	/*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) */
	/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

	/*  LWORK   (input) INTEGER */
	/*          The dimension of the array WORK. LWORK >= max(1,N). */
	/*          For optimum performance LWORK >= N*NB, where NB is the */
	/*          optimal blocksize. */

	/*          If LWORK = -1, then a workspace query is assumed; the routine */
	/*          only calculates the optimal size of the WORK array, returns */
	/*          this value as the first entry of the WORK array, and no error */
	/*          message related to LWORK is issued by XERBLA. */

	/*  INFO    (output) INTEGER */
	/*          = 0:  successful exit */
	/*          < 0:  if INFO = -i, the i-th argument has an illegal value */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Test the input arguments */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--tau;
	--work;

	/* Function Body */
	*info = 0;
	nb = ilaenv_ (&c__1, "DORGQR", " ", m, n, k, &c_n1, (ftnlen) 6, (ftnlen) 1);
	lwkopt = max (1, *n) * nb;
	work[1] = (doublereal) lwkopt;
	lquery = *lwork == -1;
	if (*m < 0)
	{
		*info = -1;
	}
	else if (*n < 0 || *n > *m)
	{
		*info = -2;
	}
	else if (*k < 0 || *k > *n)
	{
		*info = -3;
	}
	else if (*lda < max (1, *m))
	{
		*info = -5;
	}
	else if (*lwork < max (1, *n) && !lquery)
	{
		*info = -8;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DORGQR", &i__1, (ftnlen) 6);
		return 0;
	}
	else if (lquery)
	{
		return 0;
	}

	/*     Quick return if possible */

	if (*n <= 0)
	{
		work[1] = 1.;
		return 0;
	}

	nbmin = 2;
	nx = 0;
	iws = *n;
	if (nb > 1 && nb < *k)
	{

		/*        Determine when to cross over from blocked to unblocked code. */

		/* Computing MAX */
		i__1 = 0, i__2 = ilaenv_ (&c__3, "DORGQR", " ", m, n, k, &c_n1, (ftnlen) 6, (ftnlen) 1);
		nx = max (i__1, i__2);
		if (nx < *k)
		{

			/*           Determine if workspace is large enough for blocked code. */

			ldwork = *n;
			iws = ldwork * nb;
			if (*lwork < iws)
			{

				/*              Not enough workspace to use optimal NB:  reduce NB and */
				/*              determine the minimum value of NB. */

				nb = *lwork / ldwork;
				/* Computing MAX */
				i__1 = 2, i__2 = ilaenv_ (&c__2, "DORGQR", " ", m, n, k, &c_n1, (ftnlen) 6, (ftnlen) 1);
				nbmin = max (i__1, i__2);
			}
		}
	}

	if (nb >= nbmin && nb < *k && nx < *k)
	{

		/*        Use blocked code after the last block. */
		/*        The first kk columns are handled by the block method. */

		ki = (*k - nx - 1) / nb * nb;
		/* Computing MIN */
		i__1 = *k, i__2 = ki + nb;
		kk = min (i__1, i__2);

		/*        Set A(1:kk,kk+1:n) to zero. */

		i__1 = *n;
		for (j = kk + 1; j <= i__1; ++j)
		{
			i__2 = kk;
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				a[i__ + j * a_dim1] = 0.;
				/* L10: */
			}
			/* L20: */
		}
	}
	else
	{
		kk = 0;
	}

	/*     Use unblocked code for the last or only block. */

	if (kk < *n)
	{
		i__1 = *m - kk;
		i__2 = *n - kk;
		i__3 = *k - kk;
		dorg2r_ (&i__1, &i__2, &i__3, &a[kk + 1 + (kk + 1) * a_dim1], lda, &tau[kk + 1], &work[1], &iinfo);
	}

	if (kk > 0)
	{

		/*        Use blocked code */

		i__1 = -nb;
		for (i__ = ki + 1; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1)
		{
			/* Computing MIN */
			i__2 = nb, i__3 = *k - i__ + 1;
			ib = min (i__2, i__3);
			if (i__ + ib <= *n)
			{

				/*              Form the triangular factor of the block reflector */
				/*              H = H(i) H(i+1) . . . H(i+ib-1) */

				i__2 = *m - i__ + 1;
				dlarft_ ("Forward", "Columnwise", &i__2, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1], &ldwork, (ftnlen) 7, (ftnlen) 10);

				/*              Apply H to A(i:m,i+ib:n) from the left */

				i__2 = *m - i__ + 1;
				i__3 = *n - i__ - ib + 1;
				dlarfb_ ("Left", "No transpose", "Forward", "Columnwise", &i__2, &i__3, &ib, &a[i__ + i__ * a_dim1], lda, &work[1], &ldwork,
							&a[i__ + (i__ + ib) * a_dim1], lda, &work[ib + 1], &ldwork, (ftnlen) 4, (ftnlen) 12, (ftnlen) 7, (ftnlen) 10);
			}

			/*           Apply H to rows i:m of current block */

			i__2 = *m - i__ + 1;
			dorg2r_ (&i__2, &ib, &ib, &a[i__ + i__ * a_dim1], lda, &tau[i__], &work[1], &iinfo);

			/*           Set rows 1:i-1 of current block to zero */

			i__2 = i__ + ib - 1;
			for (j = i__; j <= i__2; ++j)
			{
				i__3 = i__ - 1;
				for (l = 1; l <= i__3; ++l)
				{
					a[l + j * a_dim1] = 0.;
					/* L30: */
				}
				/* L40: */
			}
			/* L50: */
		}
	}

	work[1] = (doublereal) iws;
	return 0;

	/*     End of DORGQR */

}	/* dorgqr_ */

static int
dorm2r_ (char *side, char *trans, integer * m, integer * n,
			integer * k, doublereal * a, integer * lda, doublereal * tau, doublereal *
			c__, integer * ldc, doublereal * work, integer * info, ftnlen side_len, ftnlen trans_len)
{
	/* System generated locals */
	integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;

	/* Local variables */
	static integer i__, i1, i2, i3, ic, jc, mi, ni, nq;
	static doublereal aii;
	static logical left;
	static logical notran;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     February 29, 1992 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DORM2R overwrites the general real m by n matrix C with */

	/*        Q * C  if SIDE = 'L' and TRANS = 'N', or */

	/*        Q'* C  if SIDE = 'L' and TRANS = 'T', or */

	/*        C * Q  if SIDE = 'R' and TRANS = 'N', or */

	/*        C * Q' if SIDE = 'R' and TRANS = 'T', */

	/*  where Q is a real orthogonal matrix defined as the product of k */
	/*  elementary reflectors */

	/*        Q = H(1) H(2) . . . H(k) */

	/*  as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n */
	/*  if SIDE = 'R'. */

	/*  Arguments */
	/*  ========= */

	/*  SIDE    (input) CHARACTER*1 */
	/*          = 'L': apply Q or Q' from the Left */
	/*          = 'R': apply Q or Q' from the Right */

	/*  TRANS   (input) CHARACTER*1 */
	/*          = 'N': apply Q  (No transpose) */
	/*          = 'T': apply Q' (Transpose) */

	/*  M       (input) INTEGER */
	/*          The number of rows of the matrix C. M >= 0. */

	/*  N       (input) INTEGER */
	/*          The number of columns of the matrix C. N >= 0. */

	/*  K       (input) INTEGER */
	/*          The number of elementary reflectors whose product defines */
	/*          the matrix Q. */
	/*          If SIDE = 'L', M >= K >= 0; */
	/*          if SIDE = 'R', N >= K >= 0. */

	/*  A       (input) DOUBLE PRECISION array, dimension (LDA,K) */
	/*          The i-th column must contain the vector which defines the */
	/*          elementary reflector H(i), for i = 1,2,...,k, as returned by */
	/*          DGEQRF in the first k columns of its array argument A. */
	/*          A is modified by the routine but restored on exit. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A. */
	/*          If SIDE = 'L', LDA >= max(1,M); */
	/*          if SIDE = 'R', LDA >= max(1,N). */

	/*  TAU     (input) DOUBLE PRECISION array, dimension (K) */
	/*          TAU(i) must contain the scalar factor of the elementary */
	/*          reflector H(i), as returned by DGEQRF. */

	/*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
	/*          On entry, the m by n matrix C. */
	/*          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q. */

	/*  LDC     (input) INTEGER */
	/*          The leading dimension of the array C. LDC >= max(1,M). */

	/*  WORK    (workspace) DOUBLE PRECISION array, dimension */
	/*                                   (N) if SIDE = 'L', */
	/*                                   (M) if SIDE = 'R' */

	/*  INFO    (output) INTEGER */
	/*          = 0: successful exit */
	/*          < 0: if INFO = -i, the i-th argument had an illegal value */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Test the input arguments */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	--tau;
	c_dim1 = *ldc;
	c_offset = 1 + c_dim1;
	c__ -= c_offset;
	--work;

	/* Function Body */
	*info = 0;
	left = lsame_ (side, "L", (ftnlen) 1, (ftnlen) 1);
	notran = lsame_ (trans, "N", (ftnlen) 1, (ftnlen) 1);

	/*     NQ is the order of Q */

	if (left)
	{
		nq = *m;
	}
	else
	{
		nq = *n;
	}
	if (!left && !lsame_ (side, "R", (ftnlen) 1, (ftnlen) 1))
	{
		*info = -1;
	}
	else if (!notran && !lsame_ (trans, "T", (ftnlen) 1, (ftnlen) 1))
	{
		*info = -2;
	}
	else if (*m < 0)
	{
		*info = -3;
	}
	else if (*n < 0)
	{
		*info = -4;
	}
	else if (*k < 0 || *k > nq)
	{
		*info = -5;
	}
	else if (*lda < max (1, nq))
	{
		*info = -7;
	}
	else if (*ldc < max (1, *m))
	{
		*info = -10;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DORM2R", &i__1, (ftnlen) 6);
		return 0;
	}

	/*     Quick return if possible */

	if (*m == 0 || *n == 0 || *k == 0)
	{
		return 0;
	}

	if (left && !notran || !left && notran)
	{
		i1 = 1;
		i2 = *k;
		i3 = 1;
	}
	else
	{
		i1 = *k;
		i2 = 1;
		i3 = -1;
	}

	if (left)
	{
		ni = *n;
		jc = 1;
	}
	else
	{
		mi = *m;
		ic = 1;
	}

	i__1 = i2;
	i__2 = i3;
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
	{
		if (left)
		{

			/*           H(i) is applied to C(i:m,1:n) */

			mi = *m - i__ + 1;
			ic = i__;
		}
		else
		{

			/*           H(i) is applied to C(1:m,i:n) */

			ni = *n - i__ + 1;
			jc = i__;
		}

		/*        Apply H(i) */

		aii = a[i__ + i__ * a_dim1];
		a[i__ + i__ * a_dim1] = 1.;
		dlarf_ (side, &mi, &ni, &a[i__ + i__ * a_dim1], &c__1, &tau[i__], &c__[ic + jc * c_dim1], ldc, &work[1], (ftnlen) 1);
		a[i__ + i__ * a_dim1] = aii;
		/* L10: */
	}
	return 0;

	/*     End of DORM2R */

}	/* dorm2r_ */

static int
drscl_ (integer * n, doublereal * sa, doublereal * sx, integer * incx)
{
	static doublereal mul, cden;
	static logical done;
	static doublereal cnum, cden1, cnum1;
	static doublereal bignum, smlnum;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     September 30, 1994 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DRSCL multiplies an n-element real vector x by the real scalar 1/a. */
	/*  This is done without overflow or underflow as long as */
	/*  the final result x/a does not overflow or underflow. */

	/*  Arguments */
	/*  ========= */

	/*  N       (input) INTEGER */
	/*          The number of components of the vector x. */

	/*  SA      (input) DOUBLE PRECISION */
	/*          The scalar a which is used to divide each component of x. */
	/*          SA must be >= 0, or the subroutine will divide by zero. */

	/*  SX      (input/output) DOUBLE PRECISION array, dimension */
	/*                         (1+(N-1)*abs(INCX)) */
	/*          The n-element vector x. */

	/*  INCX    (input) INTEGER */
	/*          The increment between successive values of the vector SX. */
	/*          > 0:  SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i),     1< i<= n */

	/* ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Quick return if possible */

	/* Parameter adjustments */
	--sx;

	/* Function Body */
	if (*n <= 0)
	{
		return 0;
	}

	/*     Get machine parameters */

	smlnum = dlamch_ ("S", (ftnlen) 1);
	bignum = 1. / smlnum;
	dlabad_ (&smlnum, &bignum);

	/*     Initialize the denominator to SA and the numerator to 1. */

	cden = *sa;
	cnum = 1.;

 L10:
	cden1 = cden * smlnum;
	cnum1 = cnum / bignum;
	if (abs (cden1) > abs (cnum) && cnum != 0.)
	{

		/*        Pre-multiply X by SMLNUM if CDEN is large compared to CNUM. */

		mul = smlnum;
		done = FALSE_;
		cden = cden1;
	}
	else if (abs (cnum1) > abs (cden))
	{

		/*        Pre-multiply X by BIGNUM if CDEN is small compared to CNUM. */

		mul = bignum;
		done = FALSE_;
		cnum = cnum1;
	}
	else
	{

		/*        Multiply X by CNUM / CDEN and return. */

		mul = cnum / cden;
		done = TRUE_;
	}

	/*     Scale the vector X by MUL */

	dscal_ (n, &mul, &sx[1], incx);

	if (!done)
	{
		goto L10;
	}

	return 0;

	/*     End of DRSCL */

}	/* drscl_ */

static int
dtrevc_ (char *side, char *howmny, logical * select,
			integer * n, doublereal * t, integer * ldt, doublereal * vl, integer *
			ldvl, doublereal * vr, integer * ldvr, integer * mm, integer * m, doublereal * work, integer * info, ftnlen side_len, ftnlen howmny_len)
{
	/* System generated locals */
	integer t_dim1, t_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, i__2, i__3;
	doublereal d__1, d__2, d__3, d__4;

	/* Builtin functions */
	double sqrt (doublereal);

	/* Local variables */
	static integer i__, j, k;
	static doublereal x[4] /* was [2][2] */ ;
	static integer j1, j2, n2, ii, ki, ip, is;
	static doublereal wi, wr, rec, ulp, beta, emax;
	static logical pair;
	static logical allv;
	static integer ierr;
	static doublereal unfl, ovfl, smin;
	static logical over;
	static doublereal vmax;
	static integer jnxt;
	static doublereal scale;
	static doublereal remax;
	static logical leftv, bothv;
	static doublereal vcrit;
	static logical somev;
	static doublereal xnorm;
	static doublereal bignum;
	static logical rightv;
	static doublereal smlnum;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     June 30, 1999 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DTREVC computes some or all of the right and/or left eigenvectors of */
	/*  a real upper quasi-triangular matrix T. */

	/*  The right eigenvector x and the left eigenvector y of T corresponding */
	/*  to an eigenvalue w are defined by: */

	/*               T*x = w*x,     y'*T = w*y' */

	/*  where y' denotes the conjugate transpose of the vector y. */

	/*  If all eigenvectors are requested, the routine may either return the */
	/*  matrices X and/or Y of right or left eigenvectors of T, or the */
	/*  products Q*X and/or Q*Y, where Q is an input orthogonal */
	/*  matrix. If T was obtained from the real-Schur factorization of an */
	/*  original matrix A = Q*T*Q', then Q*X and Q*Y are the matrices of */
	/*  right or left eigenvectors of A. */

	/*  T must be in Schur canonical form (as returned by DHSEQR), that is, */
	/*  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each */
	/*  2-by-2 diagonal block has its diagonal elements equal and its */
	/*  off-diagonal elements of opposite sign.  Corresponding to each 2-by-2 */
	/*  diagonal block is a complex conjugate pair of eigenvalues and */
	/*  eigenvectors; only one eigenvector of the pair is computed, namely */
	/*  the one corresponding to the eigenvalue with positive imaginary part. */

	/*  Arguments */
	/*  ========= */

	/*  SIDE    (input) CHARACTER*1 */
	/*          = 'R':  compute right eigenvectors only; */
	/*          = 'L':  compute left eigenvectors only; */
	/*          = 'B':  compute both right and left eigenvectors. */

	/*  HOWMNY  (input) CHARACTER*1 */
	/*          = 'A':  compute all right and/or left eigenvectors; */
	/*          = 'B':  compute all right and/or left eigenvectors, */
	/*                  and backtransform them using the input matrices */
	/*                  supplied in VR and/or VL; */
	/*          = 'S':  compute selected right and/or left eigenvectors, */
	/*                  specified by the logical array SELECT. */

	/*  SELECT  (input/output) LOGICAL array, dimension (N) */
	/*          If HOWMNY = 'S', SELECT specifies the eigenvectors to be */
	/*          computed. */
	/*          If HOWMNY = 'A' or 'B', SELECT is not referenced. */
	/*          To select the real eigenvector corresponding to a real */
	/*          eigenvalue w(j), SELECT(j) must be set to .TRUE..  To select */
	/*          the complex eigenvector corresponding to a complex conjugate */
	/*          pair w(j) and w(j+1), either SELECT(j) or SELECT(j+1) must be */
	/*          set to .TRUE.; then on exit SELECT(j) is .TRUE. and */
	/*          SELECT(j+1) is .FALSE.. */

	/*  N       (input) INTEGER */
	/*          The order of the matrix T. N >= 0. */

	/*  T       (input) DOUBLE PRECISION array, dimension (LDT,N) */
	/*          The upper quasi-triangular matrix T in Schur canonical form. */

	/*  LDT     (input) INTEGER */
	/*          The leading dimension of the array T. LDT >= max(1,N). */

	/*  VL      (input/output) DOUBLE PRECISION array, dimension (LDVL,MM) */
	/*          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must */
	/*          contain an N-by-N matrix Q (usually the orthogonal matrix Q */
	/*          of Schur vectors returned by DHSEQR). */
	/*          On exit, if SIDE = 'L' or 'B', VL contains: */
	/*          if HOWMNY = 'A', the matrix Y of left eigenvectors of T; */
	/*                           VL has the same quasi-lower triangular form */
	/*                           as T'. If T(i,i) is a real eigenvalue, then */
	/*                           the i-th column VL(i) of VL  is its */
	/*                           corresponding eigenvector. If T(i:i+1,i:i+1) */
	/*                           is a 2-by-2 block whose eigenvalues are */
	/*                           complex-conjugate eigenvalues of T, then */
	/*                           VL(i)+sqrt(-1)*VL(i+1) is the complex */
	/*                           eigenvector corresponding to the eigenvalue */
	/*                           with positive real part. */
	/*          if HOWMNY = 'B', the matrix Q*Y; */
	/*          if HOWMNY = 'S', the left eigenvectors of T specified by */
	/*                           SELECT, stored consecutively in the columns */
	/*                           of VL, in the same order as their */
	/*                           eigenvalues. */
	/*          A complex eigenvector corresponding to a complex eigenvalue */
	/*          is stored in two consecutive columns, the first holding the */
	/*          real part, and the second the imaginary part. */
	/*          If SIDE = 'R', VL is not referenced. */

	/*  LDVL    (input) INTEGER */
	/*          The leading dimension of the array VL.  LDVL >= max(1,N) if */
	/*          SIDE = 'L' or 'B'; LDVL >= 1 otherwise. */

	/*  VR      (input/output) DOUBLE PRECISION array, dimension (LDVR,MM) */
	/*          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must */
	/*          contain an N-by-N matrix Q (usually the orthogonal matrix Q */
	/*          of Schur vectors returned by DHSEQR). */
	/*          On exit, if SIDE = 'R' or 'B', VR contains: */
	/*          if HOWMNY = 'A', the matrix X of right eigenvectors of T; */
	/*                           VR has the same quasi-upper triangular form */
	/*                           as T. If T(i,i) is a real eigenvalue, then */
	/*                           the i-th column VR(i) of VR  is its */
	/*                           corresponding eigenvector. If T(i:i+1,i:i+1) */
	/*                           is a 2-by-2 block whose eigenvalues are */
	/*                           complex-conjugate eigenvalues of T, then */
	/*                           VR(i)+sqrt(-1)*VR(i+1) is the complex */
	/*                           eigenvector corresponding to the eigenvalue */
	/*                           with positive real part. */
	/*          if HOWMNY = 'B', the matrix Q*X; */
	/*          if HOWMNY = 'S', the right eigenvectors of T specified by */
	/*                           SELECT, stored consecutively in the columns */
	/*                           of VR, in the same order as their */
	/*                           eigenvalues. */
	/*          A complex eigenvector corresponding to a complex eigenvalue */
	/*          is stored in two consecutive columns, the first holding the */
	/*          real part and the second the imaginary part. */
	/*          If SIDE = 'L', VR is not referenced. */

	/*  LDVR    (input) INTEGER */
	/*          The leading dimension of the array VR.  LDVR >= max(1,N) if */
	/*          SIDE = 'R' or 'B'; LDVR >= 1 otherwise. */

	/*  MM      (input) INTEGER */
	/*          The number of columns in the arrays VL and/or VR. MM >= M. */

	/*  M       (output) INTEGER */
	/*          The number of columns in the arrays VL and/or VR actually */
	/*          used to store the eigenvectors. */
	/*          If HOWMNY = 'A' or 'B', M is set to N. */
	/*          Each selected real eigenvector occupies one column and each */
	/*          selected complex eigenvector occupies two columns. */

	/*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N) */

	/*  INFO    (output) INTEGER */
	/*          = 0:  successful exit */
	/*          < 0:  if INFO = -i, the i-th argument had an illegal value */

	/*  Further Details */
	/*  =============== */

	/*  The algorithm used in this program is basically backward (forward) */
	/*  substitution, with scaling to make the the code robust against */
	/*  possible overflow. */

	/*  Each eigenvector is normalized so that the element of largest */
	/*  magnitude has magnitude 1; here the magnitude of a complex number */
	/*  (x,y) is taken to be |x| + |y|. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Local Arrays .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Decode and test the input parameters */

	/* Parameter adjustments */
	--select;
	t_dim1 = *ldt;
	t_offset = 1 + t_dim1;
	t -= t_offset;
	vl_dim1 = *ldvl;
	vl_offset = 1 + vl_dim1;
	vl -= vl_offset;
	vr_dim1 = *ldvr;
	vr_offset = 1 + vr_dim1;
	vr -= vr_offset;
	--work;

	/* Function Body */
	bothv = lsame_ (side, "B", (ftnlen) 1, (ftnlen) 1);
	rightv = lsame_ (side, "R", (ftnlen) 1, (ftnlen) 1) || bothv;
	leftv = lsame_ (side, "L", (ftnlen) 1, (ftnlen) 1) || bothv;

	allv = lsame_ (howmny, "A", (ftnlen) 1, (ftnlen) 1);
	over = lsame_ (howmny, "B", (ftnlen) 1, (ftnlen) 1);
	somev = lsame_ (howmny, "S", (ftnlen) 1, (ftnlen) 1);

	*info = 0;
	if (!rightv && !leftv)
	{
		*info = -1;
	}
	else if (!allv && !over && !somev)
	{
		*info = -2;
	}
	else if (*n < 0)
	{
		*info = -4;
	}
	else if (*ldt < max (1, *n))
	{
		*info = -6;
	}
	else if (*ldvl < 1 || leftv && *ldvl < *n)
	{
		*info = -8;
	}
	else if (*ldvr < 1 || rightv && *ldvr < *n)
	{
		*info = -10;
	}
	else
	{

		/*        Set M to the number of columns required to store the selected */
		/*        eigenvectors, standardize the array SELECT if necessary, and */
		/*        test MM. */

		if (somev)
		{
			*m = 0;
			pair = FALSE_;
			i__1 = *n;
			for (j = 1; j <= i__1; ++j)
			{
				if (pair)
				{
					pair = FALSE_;
					select[j] = FALSE_;
				}
				else
				{
					if (j < *n)
					{
						if (t[j + 1 + j * t_dim1] == 0.)
						{
							if (select[j])
							{
								++(*m);
							}
						}
						else
						{
							pair = TRUE_;
							if (select[j] || select[j + 1])
							{
								select[j] = TRUE_;
								*m += 2;
							}
						}
					}
					else
					{
						if (select[*n])
						{
							++(*m);
						}
					}
				}
				/* L10: */
			}
		}
		else
		{
			*m = *n;
		}

		if (*mm < *m)
		{
			*info = -11;
		}
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DTREVC", &i__1, (ftnlen) 6);
		return 0;
	}

	/*     Quick return if possible. */

	if (*n == 0)
	{
		return 0;
	}

	/*     Set the constants to control overflow. */

	unfl = dlamch_ ("Safe minimum", (ftnlen) 12);
	ovfl = 1. / unfl;
	dlabad_ (&unfl, &ovfl);
	ulp = dlamch_ ("Precision", (ftnlen) 9);
	smlnum = unfl * (*n / ulp);
	bignum = (1. - ulp) / smlnum;

	/*     Compute 1-norm of each column of strictly upper triangular */
	/*     part of T to control overflow in triangular solver. */

	work[1] = 0.;
	i__1 = *n;
	for (j = 2; j <= i__1; ++j)
	{
		work[j] = 0.;
		i__2 = j - 1;
		for (i__ = 1; i__ <= i__2; ++i__)
		{
			work[j] += (d__1 = t[i__ + j * t_dim1], abs (d__1));
			/* L20: */
		}
		/* L30: */
	}

	/*     Index IP is used to specify the real or complex eigenvalue: */
	/*       IP = 0, real eigenvalue, */
	/*            1, first of conjugate complex pair: (wr,wi) */
	/*           -1, second of conjugate complex pair: (wr,wi) */

	n2 = *n << 1;

	if (rightv)
	{

		/*        Compute right eigenvectors. */

		ip = 0;
		is = *m;
		for (ki = *n; ki >= 1; --ki)
		{

			if (ip == 1)
			{
				goto L130;
			}
			if (ki == 1)
			{
				goto L40;
			}
			if (t[ki + (ki - 1) * t_dim1] == 0.)
			{
				goto L40;
			}
			ip = -1;

		 L40:
			if (somev)
			{
				if (ip == 0)
				{
					if (!select[ki])
					{
						goto L130;
					}
				}
				else
				{
					if (!select[ki - 1])
					{
						goto L130;
					}
				}
			}

			/*           Compute the KI-th eigenvalue (WR,WI). */

			wr = t[ki + ki * t_dim1];
			wi = 0.;
			if (ip != 0)
			{
				wi = sqrt ((d__1 = t[ki + (ki - 1) * t_dim1], abs (d__1))) * sqrt ((d__2 = t[ki - 1 + ki * t_dim1], abs (d__2)));
			}
			/* Computing MAX */
			d__1 = ulp * (abs (wr) + abs (wi));
			smin = max (d__1, smlnum);

			if (ip == 0)
			{

				/*              Real right eigenvector */

				work[ki + *n] = 1.;

				/*              Form right-hand side */

				i__1 = ki - 1;
				for (k = 1; k <= i__1; ++k)
				{
					work[k + *n] = -t[k + ki * t_dim1];
					/* L50: */
				}

				/*              Solve the upper quasi-triangular system: */
				/*                 (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK. */

				jnxt = ki - 1;
				for (j = ki - 1; j >= 1; --j)
				{
					if (j > jnxt)
					{
						goto L60;
					}
					j1 = j;
					j2 = j;
					jnxt = j - 1;
					if (j > 1)
					{
						if (t[j + (j - 1) * t_dim1] != 0.)
						{
							j1 = j - 1;
							jnxt = j - 2;
						}
					}

					if (j1 == j2)
					{

						/*                    1-by-1 diagonal block */

						dlaln2_ (&c_false, &c__1, &c__1, &smin, &c_b348, &t[j
																							 + j * t_dim1], ldt, &c_b348, &c_b348, &work[j
																																						+ *n], n, &wr, &c_b507, x, &c__2, &scale,
									&xnorm, &ierr);

						/*                    Scale X(1,1) to avoid overflow when updating */
						/*                    the right-hand side. */

						if (xnorm > 1.)
						{
							if (work[j] > bignum / xnorm)
							{
								x[0] /= xnorm;
								scale /= xnorm;
							}
						}

						/*                    Scale if necessary */

						if (scale != 1.)
						{
							dscal_ (&ki, &scale, &work[*n + 1], &c__1);
						}
						work[j + *n] = x[0];

						/*                    Update right-hand side */

						i__1 = j - 1;
						d__1 = -x[0];
						daxpy_ (&i__1, &d__1, &t[j * t_dim1 + 1], &c__1, &work[*n + 1], &c__1);

					}
					else
					{

						/*                    2-by-2 diagonal block */

						dlaln2_ (&c_false, &c__2, &c__1, &smin, &c_b348, &t[j
																							 - 1 + (j - 1) * t_dim1], ldt, &c_b348, &c_b348, &work[j - 1 + *n], n, &wr, &c_b507, x,
									&c__2, &scale, &xnorm, &ierr);

						/*                    Scale X(1,1) and X(2,1) to avoid overflow when */
						/*                    updating the right-hand side. */

						if (xnorm > 1.)
						{
							/* Computing MAX */
							d__1 = work[j - 1], d__2 = work[j];
							beta = max (d__1, d__2);
							if (beta > bignum / xnorm)
							{
								x[0] /= xnorm;
								x[1] /= xnorm;
								scale /= xnorm;
							}
						}

						/*                    Scale if necessary */

						if (scale != 1.)
						{
							dscal_ (&ki, &scale, &work[*n + 1], &c__1);
						}
						work[j - 1 + *n] = x[0];
						work[j + *n] = x[1];

						/*                    Update right-hand side */

						i__1 = j - 2;
						d__1 = -x[0];
						daxpy_ (&i__1, &d__1, &t[(j - 1) * t_dim1 + 1], &c__1, &work[*n + 1], &c__1);
						i__1 = j - 2;
						d__1 = -x[1];
						daxpy_ (&i__1, &d__1, &t[j * t_dim1 + 1], &c__1, &work[*n + 1], &c__1);
					}
				 L60:
					;
				}

				/*              Copy the vector x or Q*x to VR and normalize. */

				if (!over)
				{
					dcopy_ (&ki, &work[*n + 1], &c__1, &vr[is * vr_dim1 + 1], &c__1);

					ii = idamax_ (&ki, &vr[is * vr_dim1 + 1], &c__1);
					remax = 1. / (d__1 = vr[ii + is * vr_dim1], abs (d__1));
					dscal_ (&ki, &remax, &vr[is * vr_dim1 + 1], &c__1);

					i__1 = *n;
					for (k = ki + 1; k <= i__1; ++k)
					{
						vr[k + is * vr_dim1] = 0.;
						/* L70: */
					}
				}
				else
				{
					if (ki > 1)
					{
						i__1 = ki - 1;
						dgemv_ ("N", n, &i__1, &c_b348, &vr[vr_offset], ldvr, &work[*n + 1], &c__1, &work[ki + *n], &vr[ki * vr_dim1 + 1], &c__1, (ftnlen) 1);
					}

					ii = idamax_ (n, &vr[ki * vr_dim1 + 1], &c__1);
					remax = 1. / (d__1 = vr[ii + ki * vr_dim1], abs (d__1));
					dscal_ (n, &remax, &vr[ki * vr_dim1 + 1], &c__1);
				}

			}
			else
			{

				/*              Complex right eigenvector. */

				/*              Initial solve */
				/*                [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = 0. */
				/*                [ (T(KI,KI-1)   T(KI,KI)   )               ] */

				if ((d__1 = t[ki - 1 + ki * t_dim1], abs (d__1)) >= (d__2 = t[ki + (ki - 1) * t_dim1], abs (d__2)))
				{
					work[ki - 1 + *n] = 1.;
					work[ki + n2] = wi / t[ki - 1 + ki * t_dim1];
				}
				else
				{
					work[ki - 1 + *n] = -wi / t[ki + (ki - 1) * t_dim1];
					work[ki + n2] = 1.;
				}
				work[ki + *n] = 0.;
				work[ki - 1 + n2] = 0.;

				/*              Form right-hand side */

				i__1 = ki - 2;
				for (k = 1; k <= i__1; ++k)
				{
					work[k + *n] = -work[ki - 1 + *n] * t[k + (ki - 1) * t_dim1];
					work[k + n2] = -work[ki + n2] * t[k + ki * t_dim1];
					/* L80: */
				}

				/*              Solve upper quasi-triangular system: */
				/*              (T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2) */

				jnxt = ki - 2;
				for (j = ki - 2; j >= 1; --j)
				{
					if (j > jnxt)
					{
						goto L90;
					}
					j1 = j;
					j2 = j;
					jnxt = j - 1;
					if (j > 1)
					{
						if (t[j + (j - 1) * t_dim1] != 0.)
						{
							j1 = j - 1;
							jnxt = j - 2;
						}
					}

					if (j1 == j2)
					{

						/*                    1-by-1 diagonal block */

						dlaln2_ (&c_false, &c__1, &c__2, &smin, &c_b348, &t[j
																							 + j * t_dim1], ldt, &c_b348, &c_b348, &work[j
																																						+ *n], n, &wr, &wi, x, &c__2, &scale, &xnorm,
									&ierr);

						/*                    Scale X(1,1) and X(1,2) to avoid overflow when */
						/*                    updating the right-hand side. */

						if (xnorm > 1.)
						{
							if (work[j] > bignum / xnorm)
							{
								x[0] /= xnorm;
								x[2] /= xnorm;
								scale /= xnorm;
							}
						}

						/*                    Scale if necessary */

						if (scale != 1.)
						{
							dscal_ (&ki, &scale, &work[*n + 1], &c__1);
							dscal_ (&ki, &scale, &work[n2 + 1], &c__1);
						}
						work[j + *n] = x[0];
						work[j + n2] = x[2];

						/*                    Update the right-hand side */

						i__1 = j - 1;
						d__1 = -x[0];
						daxpy_ (&i__1, &d__1, &t[j * t_dim1 + 1], &c__1, &work[*n + 1], &c__1);
						i__1 = j - 1;
						d__1 = -x[2];
						daxpy_ (&i__1, &d__1, &t[j * t_dim1 + 1], &c__1, &work[n2 + 1], &c__1);

					}
					else
					{

						/*                    2-by-2 diagonal block */

						dlaln2_ (&c_false, &c__2, &c__2, &smin, &c_b348, &t[j
																							 - 1 + (j - 1) * t_dim1], ldt, &c_b348, &c_b348, &work[j - 1 + *n], n, &wr, &wi, x, &c__2,
									&scale, &xnorm, &ierr);

						/*                    Scale X to avoid overflow when updating */
						/*                    the right-hand side. */

						if (xnorm > 1.)
						{
							/* Computing MAX */
							d__1 = work[j - 1], d__2 = work[j];
							beta = max (d__1, d__2);
							if (beta > bignum / xnorm)
							{
								rec = 1. / xnorm;
								x[0] *= rec;
								x[2] *= rec;
								x[1] *= rec;
								x[3] *= rec;
								scale *= rec;
							}
						}

						/*                    Scale if necessary */

						if (scale != 1.)
						{
							dscal_ (&ki, &scale, &work[*n + 1], &c__1);
							dscal_ (&ki, &scale, &work[n2 + 1], &c__1);
						}
						work[j - 1 + *n] = x[0];
						work[j + *n] = x[1];
						work[j - 1 + n2] = x[2];
						work[j + n2] = x[3];

						/*                    Update the right-hand side */

						i__1 = j - 2;
						d__1 = -x[0];
						daxpy_ (&i__1, &d__1, &t[(j - 1) * t_dim1 + 1], &c__1, &work[*n + 1], &c__1);
						i__1 = j - 2;
						d__1 = -x[1];
						daxpy_ (&i__1, &d__1, &t[j * t_dim1 + 1], &c__1, &work[*n + 1], &c__1);
						i__1 = j - 2;
						d__1 = -x[2];
						daxpy_ (&i__1, &d__1, &t[(j - 1) * t_dim1 + 1], &c__1, &work[n2 + 1], &c__1);
						i__1 = j - 2;
						d__1 = -x[3];
						daxpy_ (&i__1, &d__1, &t[j * t_dim1 + 1], &c__1, &work[n2 + 1], &c__1);
					}
				 L90:
					;
				}

				/*              Copy the vector x or Q*x to VR and normalize. */

				if (!over)
				{
					dcopy_ (&ki, &work[*n + 1], &c__1, &vr[(is - 1) * vr_dim1 + 1], &c__1);
					dcopy_ (&ki, &work[n2 + 1], &c__1, &vr[is * vr_dim1 + 1], &c__1);

					emax = 0.;
					i__1 = ki;
					for (k = 1; k <= i__1; ++k)
					{
						/* Computing MAX */
						d__3 = emax, d__4 = (d__1 = vr[k + (is - 1) * vr_dim1], abs (d__1)) + (d__2 = vr[k + is * vr_dim1], abs (d__2));
						emax = max (d__3, d__4);
						/* L100: */
					}

					remax = 1. / emax;
					dscal_ (&ki, &remax, &vr[(is - 1) * vr_dim1 + 1], &c__1);
					dscal_ (&ki, &remax, &vr[is * vr_dim1 + 1], &c__1);

					i__1 = *n;
					for (k = ki + 1; k <= i__1; ++k)
					{
						vr[k + (is - 1) * vr_dim1] = 0.;
						vr[k + is * vr_dim1] = 0.;
						/* L110: */
					}

				}
				else
				{

					if (ki > 2)
					{
						i__1 = ki - 2;
						dgemv_ ("N", n, &i__1, &c_b348, &vr[vr_offset], ldvr, &work[*n + 1], &c__1, &work[ki - 1 + *n], &vr[(ki - 1) * vr_dim1 + 1], &c__1,
								  (ftnlen) 1);
						i__1 = ki - 2;
						dgemv_ ("N", n, &i__1, &c_b348, &vr[vr_offset], ldvr, &work[n2 + 1], &c__1, &work[ki + n2], &vr[ki * vr_dim1 + 1], &c__1, (ftnlen) 1);
					}
					else
					{
						dscal_ (n, &work[ki - 1 + *n], &vr[(ki - 1) * vr_dim1 + 1], &c__1);
						dscal_ (n, &work[ki + n2], &vr[ki * vr_dim1 + 1], &c__1);
					}

					emax = 0.;
					i__1 = *n;
					for (k = 1; k <= i__1; ++k)
					{
						/* Computing MAX */
						d__3 = emax, d__4 = (d__1 = vr[k + (ki - 1) * vr_dim1], abs (d__1)) + (d__2 = vr[k + ki * vr_dim1], abs (d__2));
						emax = max (d__3, d__4);
						/* L120: */
					}
					remax = 1. / emax;
					dscal_ (n, &remax, &vr[(ki - 1) * vr_dim1 + 1], &c__1);
					dscal_ (n, &remax, &vr[ki * vr_dim1 + 1], &c__1);
				}
			}

			--is;
			if (ip != 0)
			{
				--is;
			}
		 L130:
			if (ip == 1)
			{
				ip = 0;
			}
			if (ip == -1)
			{
				ip = 1;
			}
			/* L140: */
		}
	}

	if (leftv)
	{

		/*        Compute left eigenvectors. */

		ip = 0;
		is = 1;
		i__1 = *n;
		for (ki = 1; ki <= i__1; ++ki)
		{

			if (ip == -1)
			{
				goto L250;
			}
			if (ki == *n)
			{
				goto L150;
			}
			if (t[ki + 1 + ki * t_dim1] == 0.)
			{
				goto L150;
			}
			ip = 1;

		 L150:
			if (somev)
			{
				if (!select[ki])
				{
					goto L250;
				}
			}

			/*           Compute the KI-th eigenvalue (WR,WI). */

			wr = t[ki + ki * t_dim1];
			wi = 0.;
			if (ip != 0)
			{
				wi = sqrt ((d__1 = t[ki + (ki + 1) * t_dim1], abs (d__1))) * sqrt ((d__2 = t[ki + 1 + ki * t_dim1], abs (d__2)));
			}
			/* Computing MAX */
			d__1 = ulp * (abs (wr) + abs (wi));
			smin = max (d__1, smlnum);

			if (ip == 0)
			{

				/*              Real left eigenvector. */

				work[ki + *n] = 1.;

				/*              Form right-hand side */

				i__2 = *n;
				for (k = ki + 1; k <= i__2; ++k)
				{
					work[k + *n] = -t[ki + k * t_dim1];
					/* L160: */
				}

				/*              Solve the quasi-triangular system: */
				/*                 (T(KI+1:N,KI+1:N) - WR)'*X = SCALE*WORK */

				vmax = 1.;
				vcrit = bignum;

				jnxt = ki + 1;
				i__2 = *n;
				for (j = ki + 1; j <= i__2; ++j)
				{
					if (j < jnxt)
					{
						goto L170;
					}
					j1 = j;
					j2 = j;
					jnxt = j + 1;
					if (j < *n)
					{
						if (t[j + 1 + j * t_dim1] != 0.)
						{
							j2 = j + 1;
							jnxt = j + 2;
						}
					}

					if (j1 == j2)
					{

						/*                    1-by-1 diagonal block */

						/*                    Scale if necessary to avoid overflow when forming */
						/*                    the right-hand side. */

						if (work[j] > vcrit)
						{
							rec = 1. / vmax;
							i__3 = *n - ki + 1;
							dscal_ (&i__3, &rec, &work[ki + *n], &c__1);
							vmax = 1.;
							vcrit = bignum;
						}

						i__3 = j - ki - 1;
						work[j + *n] -= ddot_ (&i__3, &t[ki + 1 + j * t_dim1], &c__1, &work[ki + 1 + *n], &c__1);

						/*                    Solve (T(J,J)-WR)'*X = WORK */

						dlaln2_ (&c_false, &c__1, &c__1, &smin, &c_b348, &t[j
																							 + j * t_dim1], ldt, &c_b348, &c_b348, &work[j
																																						+ *n], n, &wr, &c_b507, x, &c__2, &scale,
									&xnorm, &ierr);

						/*                    Scale if necessary */

						if (scale != 1.)
						{
							i__3 = *n - ki + 1;
							dscal_ (&i__3, &scale, &work[ki + *n], &c__1);
						}
						work[j + *n] = x[0];
						/* Computing MAX */
						d__2 = (d__1 = work[j + *n], abs (d__1));
						vmax = max (d__2, vmax);
						vcrit = bignum / vmax;

					}
					else
					{

						/*                    2-by-2 diagonal block */

						/*                    Scale if necessary to avoid overflow when forming */
						/*                    the right-hand side. */

						/* Computing MAX */
						d__1 = work[j], d__2 = work[j + 1];
						beta = max (d__1, d__2);
						if (beta > vcrit)
						{
							rec = 1. / vmax;
							i__3 = *n - ki + 1;
							dscal_ (&i__3, &rec, &work[ki + *n], &c__1);
							vmax = 1.;
							vcrit = bignum;
						}

						i__3 = j - ki - 1;
						work[j + *n] -= ddot_ (&i__3, &t[ki + 1 + j * t_dim1], &c__1, &work[ki + 1 + *n], &c__1);

						i__3 = j - ki - 1;
						work[j + 1 + *n] -= ddot_ (&i__3, &t[ki + 1 + (j + 1) * t_dim1], &c__1, &work[ki + 1 + *n], &c__1);

						/*                    Solve */
						/*                      [T(J,J)-WR   T(J,J+1)     ]'* X = SCALE*( WORK1 ) */
						/*                      [T(J+1,J)    T(J+1,J+1)-WR]             ( WORK2 ) */

						dlaln2_ (&c_true, &c__2, &c__1, &smin, &c_b348, &t[j +
																							j * t_dim1], ldt, &c_b348, &c_b348, &work[j +
																																					*n], n, &wr, &c_b507, x, &c__2, &scale, &xnorm,
									&ierr);

						/*                    Scale if necessary */

						if (scale != 1.)
						{
							i__3 = *n - ki + 1;
							dscal_ (&i__3, &scale, &work[ki + *n], &c__1);
						}
						work[j + *n] = x[0];
						work[j + 1 + *n] = x[1];

						/* Computing MAX */
						d__3 = (d__1 = work[j + *n], abs (d__1)), d__4 = (d__2 = work[j + 1 + *n], abs (d__2)), d__3 = max (d__3, d__4);
						vmax = max (d__3, vmax);
						vcrit = bignum / vmax;

					}
				 L170:
					;
				}

				/*              Copy the vector x or Q*x to VL and normalize. */

				if (!over)
				{
					i__2 = *n - ki + 1;
					dcopy_ (&i__2, &work[ki + *n], &c__1, &vl[ki + is * vl_dim1], &c__1);

					i__2 = *n - ki + 1;
					ii = idamax_ (&i__2, &vl[ki + is * vl_dim1], &c__1) + ki - 1;
					remax = 1. / (d__1 = vl[ii + is * vl_dim1], abs (d__1));
					i__2 = *n - ki + 1;
					dscal_ (&i__2, &remax, &vl[ki + is * vl_dim1], &c__1);

					i__2 = ki - 1;
					for (k = 1; k <= i__2; ++k)
					{
						vl[k + is * vl_dim1] = 0.;
						/* L180: */
					}

				}
				else
				{

					if (ki < *n)
					{
						i__2 = *n - ki;
						dgemv_ ("N", n, &i__2, &c_b348, &vl[(ki + 1) * vl_dim1
																		+ 1], ldvl, &work[ki + 1 + *n], &c__1, &work[ki + *n], &vl[ki * vl_dim1 + 1], &c__1, (ftnlen) 1);
					}

					ii = idamax_ (n, &vl[ki * vl_dim1 + 1], &c__1);
					remax = 1. / (d__1 = vl[ii + ki * vl_dim1], abs (d__1));
					dscal_ (n, &remax, &vl[ki * vl_dim1 + 1], &c__1);

				}

			}
			else
			{

				/*              Complex left eigenvector. */

				/*               Initial solve: */
				/*                 ((T(KI,KI)    T(KI,KI+1) )' - (WR - I* WI))*X = 0. */
				/*                 ((T(KI+1,KI) T(KI+1,KI+1))                ) */

				if ((d__1 = t[ki + (ki + 1) * t_dim1], abs (d__1)) >= (d__2 = t[ki + 1 + ki * t_dim1], abs (d__2)))
				{
					work[ki + *n] = wi / t[ki + (ki + 1) * t_dim1];
					work[ki + 1 + n2] = 1.;
				}
				else
				{
					work[ki + *n] = 1.;
					work[ki + 1 + n2] = -wi / t[ki + 1 + ki * t_dim1];
				}
				work[ki + 1 + *n] = 0.;
				work[ki + n2] = 0.;

				/*              Form right-hand side */

				i__2 = *n;
				for (k = ki + 2; k <= i__2; ++k)
				{
					work[k + *n] = -work[ki + *n] * t[ki + k * t_dim1];
					work[k + n2] = -work[ki + 1 + n2] * t[ki + 1 + k * t_dim1];
					/* L190: */
				}

				/*              Solve complex quasi-triangular system: */
				/*              ( T(KI+2,N:KI+2,N) - (WR-i*WI) )*X = WORK1+i*WORK2 */

				vmax = 1.;
				vcrit = bignum;

				jnxt = ki + 2;
				i__2 = *n;
				for (j = ki + 2; j <= i__2; ++j)
				{
					if (j < jnxt)
					{
						goto L200;
					}
					j1 = j;
					j2 = j;
					jnxt = j + 1;
					if (j < *n)
					{
						if (t[j + 1 + j * t_dim1] != 0.)
						{
							j2 = j + 1;
							jnxt = j + 2;
						}
					}

					if (j1 == j2)
					{

						/*                    1-by-1 diagonal block */

						/*                    Scale if necessary to avoid overflow when */
						/*                    forming the right-hand side elements. */

						if (work[j] > vcrit)
						{
							rec = 1. / vmax;
							i__3 = *n - ki + 1;
							dscal_ (&i__3, &rec, &work[ki + *n], &c__1);
							i__3 = *n - ki + 1;
							dscal_ (&i__3, &rec, &work[ki + n2], &c__1);
							vmax = 1.;
							vcrit = bignum;
						}

						i__3 = j - ki - 2;
						work[j + *n] -= ddot_ (&i__3, &t[ki + 2 + j * t_dim1], &c__1, &work[ki + 2 + *n], &c__1);
						i__3 = j - ki - 2;
						work[j + n2] -= ddot_ (&i__3, &t[ki + 2 + j * t_dim1], &c__1, &work[ki + 2 + n2], &c__1);

						/*                    Solve (T(J,J)-(WR-i*WI))*(X11+i*X12)= WK+I*WK2 */

						d__1 = -wi;
						dlaln2_ (&c_false, &c__1, &c__2, &smin, &c_b348, &t[j
																							 + j * t_dim1], ldt, &c_b348, &c_b348, &work[j
																																						+ *n], n, &wr, &d__1, x, &c__2, &scale,
									&xnorm, &ierr);

						/*                    Scale if necessary */

						if (scale != 1.)
						{
							i__3 = *n - ki + 1;
							dscal_ (&i__3, &scale, &work[ki + *n], &c__1);
							i__3 = *n - ki + 1;
							dscal_ (&i__3, &scale, &work[ki + n2], &c__1);
						}
						work[j + *n] = x[0];
						work[j + n2] = x[2];
						/* Computing MAX */
						d__3 = (d__1 = work[j + *n], abs (d__1)), d__4 = (d__2 = work[j + n2], abs (d__2)), d__3 = max (d__3, d__4);
						vmax = max (d__3, vmax);
						vcrit = bignum / vmax;

					}
					else
					{

						/*                    2-by-2 diagonal block */

						/*                    Scale if necessary to avoid overflow when forming */
						/*                    the right-hand side elements. */

						/* Computing MAX */
						d__1 = work[j], d__2 = work[j + 1];
						beta = max (d__1, d__2);
						if (beta > vcrit)
						{
							rec = 1. / vmax;
							i__3 = *n - ki + 1;
							dscal_ (&i__3, &rec, &work[ki + *n], &c__1);
							i__3 = *n - ki + 1;
							dscal_ (&i__3, &rec, &work[ki + n2], &c__1);
							vmax = 1.;
							vcrit = bignum;
						}

						i__3 = j - ki - 2;
						work[j + *n] -= ddot_ (&i__3, &t[ki + 2 + j * t_dim1], &c__1, &work[ki + 2 + *n], &c__1);

						i__3 = j - ki - 2;
						work[j + n2] -= ddot_ (&i__3, &t[ki + 2 + j * t_dim1], &c__1, &work[ki + 2 + n2], &c__1);

						i__3 = j - ki - 2;
						work[j + 1 + *n] -= ddot_ (&i__3, &t[ki + 2 + (j + 1) * t_dim1], &c__1, &work[ki + 2 + *n], &c__1);

						i__3 = j - ki - 2;
						work[j + 1 + n2] -= ddot_ (&i__3, &t[ki + 2 + (j + 1) * t_dim1], &c__1, &work[ki + 2 + n2], &c__1);

						/*                    Solve 2-by-2 complex linear equation */
						/*                      ([T(j,j)   T(j,j+1)  ]'-(wr-i*wi)*I)*X = SCALE*B */
						/*                      ([T(j+1,j) T(j+1,j+1)]             ) */

						d__1 = -wi;
						dlaln2_ (&c_true, &c__2, &c__2, &smin, &c_b348, &t[j +
																							j * t_dim1], ldt, &c_b348, &c_b348, &work[j +
																																					*n], n, &wr, &d__1, x, &c__2, &scale, &xnorm,
									&ierr);

						/*                    Scale if necessary */

						if (scale != 1.)
						{
							i__3 = *n - ki + 1;
							dscal_ (&i__3, &scale, &work[ki + *n], &c__1);
							i__3 = *n - ki + 1;
							dscal_ (&i__3, &scale, &work[ki + n2], &c__1);
						}
						work[j + *n] = x[0];
						work[j + n2] = x[2];
						work[j + 1 + *n] = x[1];
						work[j + 1 + n2] = x[3];
						/* Computing MAX */
						d__1 = abs (x[0]), d__2 = abs (x[2]), d__1 = max (d__1,
																						  d__2), d__2 = abs (x[1]), d__1 = max (d__1, d__2), d__2 = abs (x[3]), d__1 =
							max (d__1, d__2);
						vmax = max (d__1, vmax);
						vcrit = bignum / vmax;

					}
				 L200:
					;
				}

				/*              Copy the vector x or Q*x to VL and normalize. */

				/* L210: */
				if (!over)
				{
					i__2 = *n - ki + 1;
					dcopy_ (&i__2, &work[ki + *n], &c__1, &vl[ki + is * vl_dim1], &c__1);
					i__2 = *n - ki + 1;
					dcopy_ (&i__2, &work[ki + n2], &c__1, &vl[ki + (is + 1) * vl_dim1], &c__1);

					emax = 0.;
					i__2 = *n;
					for (k = ki; k <= i__2; ++k)
					{
						/* Computing MAX */
						d__3 = emax, d__4 = (d__1 = vl[k + is * vl_dim1], abs (d__1)) + (d__2 = vl[k + (is + 1) * vl_dim1], abs (d__2));
						emax = max (d__3, d__4);
						/* L220: */
					}
					remax = 1. / emax;
					i__2 = *n - ki + 1;
					dscal_ (&i__2, &remax, &vl[ki + is * vl_dim1], &c__1);
					i__2 = *n - ki + 1;
					dscal_ (&i__2, &remax, &vl[ki + (is + 1) * vl_dim1], &c__1);

					i__2 = ki - 1;
					for (k = 1; k <= i__2; ++k)
					{
						vl[k + is * vl_dim1] = 0.;
						vl[k + (is + 1) * vl_dim1] = 0.;
						/* L230: */
					}
				}
				else
				{
					if (ki < *n - 1)
					{
						i__2 = *n - ki - 1;
						dgemv_ ("N", n, &i__2, &c_b348, &vl[(ki + 2) * vl_dim1
																		+ 1], ldvl, &work[ki + 2 + *n], &c__1, &work[ki + *n], &vl[ki * vl_dim1 + 1], &c__1, (ftnlen) 1);
						i__2 = *n - ki - 1;
						dgemv_ ("N", n, &i__2, &c_b348, &vl[(ki + 2) * vl_dim1
																		+ 1], ldvl, &work[ki + 2 + n2], &c__1, &work[ki + 1 + n2], &vl[(ki + 1) * vl_dim1 + 1], &c__1,
								  (ftnlen) 1);
					}
					else
					{
						dscal_ (n, &work[ki + *n], &vl[ki * vl_dim1 + 1], &c__1);
						dscal_ (n, &work[ki + 1 + n2], &vl[(ki + 1) * vl_dim1 + 1], &c__1);
					}

					emax = 0.;
					i__2 = *n;
					for (k = 1; k <= i__2; ++k)
					{
						/* Computing MAX */
						d__3 = emax, d__4 = (d__1 = vl[k + ki * vl_dim1], abs (d__1)) + (d__2 = vl[k + (ki + 1) * vl_dim1], abs (d__2));
						emax = max (d__3, d__4);
						/* L240: */
					}
					remax = 1. / emax;
					dscal_ (n, &remax, &vl[ki * vl_dim1 + 1], &c__1);
					dscal_ (n, &remax, &vl[(ki + 1) * vl_dim1 + 1], &c__1);

				}

			}

			++is;
			if (ip != 0)
			{
				++is;
			}
		 L250:
			if (ip == -1)
			{
				ip = 0;
			}
			if (ip == 1)
			{
				ip = -1;
			}

			/* L260: */
		}

	}

	return 0;

	/*     End of DTREVC */

}	/* dtrevc_ */

static int
dtrexc_ (char *compq, integer * n, doublereal * t, integer *
			ldt, doublereal * q, integer * ldq, integer * ifst, integer * ilst, doublereal * work, integer * info, ftnlen compq_len)
{
	/* System generated locals */
	integer q_dim1, q_offset, t_dim1, t_offset, i__1;

	/* Local variables */
	static integer nbf, nbl, here;
	static logical wantq;
	static integer nbnext;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     March 31, 1993 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DTREXC reorders the real Schur factorization of a real matrix */
	/*  A = Q*T*Q**T, so that the diagonal block of T with row index IFST is */
	/*  moved to row ILST. */

	/*  The real Schur form T is reordered by an orthogonal similarity */
	/*  transformation Z**T*T*Z, and optionally the matrix Q of Schur vectors */
	/*  is updated by postmultiplying it with Z. */

	/*  T must be in Schur canonical form (as returned by DHSEQR), that is, */
	/*  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each */
	/*  2-by-2 diagonal block has its diagonal elements equal and its */
	/*  off-diagonal elements of opposite sign. */

	/*  Arguments */
	/*  ========= */

	/*  COMPQ   (input) CHARACTER*1 */
	/*          = 'V':  update the matrix Q of Schur vectors; */
	/*          = 'N':  do not update Q. */

	/*  N       (input) INTEGER */
	/*          The order of the matrix T. N >= 0. */

	/*  T       (input/output) DOUBLE PRECISION array, dimension (LDT,N) */
	/*          On entry, the upper quasi-triangular matrix T, in Schur */
	/*          Schur canonical form. */
	/*          On exit, the reordered upper quasi-triangular matrix, again */
	/*          in Schur canonical form. */

	/*  LDT     (input) INTEGER */
	/*          The leading dimension of the array T. LDT >= max(1,N). */

	/*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
	/*          On entry, if COMPQ = 'V', the matrix Q of Schur vectors. */
	/*          On exit, if COMPQ = 'V', Q has been postmultiplied by the */
	/*          orthogonal transformation matrix Z which reorders T. */
	/*          If COMPQ = 'N', Q is not referenced. */

	/*  LDQ     (input) INTEGER */
	/*          The leading dimension of the array Q.  LDQ >= max(1,N). */

	/*  IFST    (input/output) INTEGER */
	/*  ILST    (input/output) INTEGER */
	/*          Specify the reordering of the diagonal blocks of T. */
	/*          The block with row index IFST is moved to row ILST, by a */
	/*          sequence of transpositions between adjacent blocks. */
	/*          On exit, if IFST pointed on entry to the second row of a */
	/*          2-by-2 block, it is changed to point to the first row; ILST */
	/*          always points to the first row of the block in its final */
	/*          position (which may differ from its input value by +1 or -1). */
	/*          1 <= IFST <= N; 1 <= ILST <= N. */

	/*  WORK    (workspace) DOUBLE PRECISION array, dimension (N) */

	/*  INFO    (output) INTEGER */
	/*          = 0:  successful exit */
	/*          < 0:  if INFO = -i, the i-th argument had an illegal value */
	/*          = 1:  two adjacent blocks were too close to swap (the problem */
	/*                is very ill-conditioned); T may have been partially */
	/*                reordered, and ILST points to the first row of the */
	/*                current position of the block being moved. */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Decode and test the input arguments. */

	/* Parameter adjustments */
	t_dim1 = *ldt;
	t_offset = 1 + t_dim1;
	t -= t_offset;
	q_dim1 = *ldq;
	q_offset = 1 + q_dim1;
	q -= q_offset;
	--work;

	/* Function Body */
	*info = 0;
	wantq = lsame_ (compq, "V", (ftnlen) 1, (ftnlen) 1);
	if (!wantq && !lsame_ (compq, "N", (ftnlen) 1, (ftnlen) 1))
	{
		*info = -1;
	}
	else if (*n < 0)
	{
		*info = -2;
	}
	else if (*ldt < max (1, *n))
	{
		*info = -4;
	}
	else if (*ldq < 1 || wantq && *ldq < max (1, *n))
	{
		*info = -6;
	}
	else if (*ifst < 1 || *ifst > *n)
	{
		*info = -7;
	}
	else if (*ilst < 1 || *ilst > *n)
	{
		*info = -8;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DTREXC", &i__1, (ftnlen) 6);
		return 0;
	}

	/*     Quick return if possible */

	if (*n <= 1)
	{
		return 0;
	}

	/*     Determine the first row of specified block */
	/*     and find out it is 1 by 1 or 2 by 2. */

	if (*ifst > 1)
	{
		if (t[*ifst + (*ifst - 1) * t_dim1] != 0.)
		{
			--(*ifst);
		}
	}
	nbf = 1;
	if (*ifst < *n)
	{
		if (t[*ifst + 1 + *ifst * t_dim1] != 0.)
		{
			nbf = 2;
		}
	}

	/*     Determine the first row of the final block */
	/*     and find out it is 1 by 1 or 2 by 2. */

	if (*ilst > 1)
	{
		if (t[*ilst + (*ilst - 1) * t_dim1] != 0.)
		{
			--(*ilst);
		}
	}
	nbl = 1;
	if (*ilst < *n)
	{
		if (t[*ilst + 1 + *ilst * t_dim1] != 0.)
		{
			nbl = 2;
		}
	}

	if (*ifst == *ilst)
	{
		return 0;
	}

	if (*ifst < *ilst)
	{

		/*        Update ILST */

		if (nbf == 2 && nbl == 1)
		{
			--(*ilst);
		}
		if (nbf == 1 && nbl == 2)
		{
			++(*ilst);
		}

		here = *ifst;

	 L10:

		/*        Swap block with next one below */

		if (nbf == 1 || nbf == 2)
		{

			/*           Current block either 1 by 1 or 2 by 2 */

			nbnext = 1;
			if (here + nbf + 1 <= *n)
			{
				if (t[here + nbf + 1 + (here + nbf) * t_dim1] != 0.)
				{
					nbnext = 2;
				}
			}
			dlaexc_ (&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &here, &nbf, &nbnext, &work[1], info);
			if (*info != 0)
			{
				*ilst = here;
				return 0;
			}
			here += nbnext;

			/*           Test if 2 by 2 block breaks into two 1 by 1 blocks */

			if (nbf == 2)
			{
				if (t[here + 1 + here * t_dim1] == 0.)
				{
					nbf = 3;
				}
			}

		}
		else
		{

			/*           Current block consists of two 1 by 1 blocks each of which */
			/*           must be swapped individually */

			nbnext = 1;
			if (here + 3 <= *n)
			{
				if (t[here + 3 + (here + 2) * t_dim1] != 0.)
				{
					nbnext = 2;
				}
			}
			i__1 = here + 1;
			dlaexc_ (&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &c__1, &nbnext, &work[1], info);
			if (*info != 0)
			{
				*ilst = here;
				return 0;
			}
			if (nbnext == 1)
			{

				/*              Swap two 1 by 1 blocks, no problems possible */

				dlaexc_ (&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &here, &c__1, &nbnext, &work[1], info);
				++here;
			}
			else
			{

				/*              Recompute NBNEXT in case 2 by 2 split */

				if (t[here + 2 + (here + 1) * t_dim1] == 0.)
				{
					nbnext = 1;
				}
				if (nbnext == 2)
				{

					/*                 2 by 2 Block did not split */

					dlaexc_ (&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &here, &c__1, &nbnext, &work[1], info);
					if (*info != 0)
					{
						*ilst = here;
						return 0;
					}
					here += 2;
				}
				else
				{

					/*                 2 by 2 Block did split */

					dlaexc_ (&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &here, &c__1, &c__1, &work[1], info);
					i__1 = here + 1;
					dlaexc_ (&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &c__1, &c__1, &work[1], info);
					here += 2;
				}
			}
		}
		if (here < *ilst)
		{
			goto L10;
		}

	}
	else
	{

		here = *ifst;
	 L20:

		/*        Swap block with next one above */

		if (nbf == 1 || nbf == 2)
		{

			/*           Current block either 1 by 1 or 2 by 2 */

			nbnext = 1;
			if (here >= 3)
			{
				if (t[here - 1 + (here - 2) * t_dim1] != 0.)
				{
					nbnext = 2;
				}
			}
			i__1 = here - nbnext;
			dlaexc_ (&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &nbnext, &nbf, &work[1], info);
			if (*info != 0)
			{
				*ilst = here;
				return 0;
			}
			here -= nbnext;

			/*           Test if 2 by 2 block breaks into two 1 by 1 blocks */

			if (nbf == 2)
			{
				if (t[here + 1 + here * t_dim1] == 0.)
				{
					nbf = 3;
				}
			}

		}
		else
		{

			/*           Current block consists of two 1 by 1 blocks each of which */
			/*           must be swapped individually */

			nbnext = 1;
			if (here >= 3)
			{
				if (t[here - 1 + (here - 2) * t_dim1] != 0.)
				{
					nbnext = 2;
				}
			}
			i__1 = here - nbnext;
			dlaexc_ (&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &nbnext, &c__1, &work[1], info);
			if (*info != 0)
			{
				*ilst = here;
				return 0;
			}
			if (nbnext == 1)
			{

				/*              Swap two 1 by 1 blocks, no problems possible */

				dlaexc_ (&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &here, &nbnext, &c__1, &work[1], info);
				--here;
			}
			else
			{

				/*              Recompute NBNEXT in case 2 by 2 split */

				if (t[here + (here - 1) * t_dim1] == 0.)
				{
					nbnext = 1;
				}
				if (nbnext == 2)
				{

					/*                 2 by 2 Block did not split */

					i__1 = here - 1;
					dlaexc_ (&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &c__2, &c__1, &work[1], info);
					if (*info != 0)
					{
						*ilst = here;
						return 0;
					}
					here += -2;
				}
				else
				{

					/*                 2 by 2 Block did split */

					dlaexc_ (&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &here, &c__1, &c__1, &work[1], info);
					i__1 = here - 1;
					dlaexc_ (&wantq, n, &t[t_offset], ldt, &q[q_offset], ldq, &i__1, &c__1, &c__1, &work[1], info);
					here += -2;
				}
			}
		}
		if (here > *ilst)
		{
			goto L20;
		}
	}
	*ilst = here;

	return 0;

	/*     End of DTREXC */

}	/* dtrexc_ */

static int
dtrsen_ (char *job, char *compq, logical * select, integer
			* n, doublereal * t, integer * ldt, doublereal * q, integer * ldq,
			doublereal * wr, doublereal * wi, integer * m, doublereal * s, doublereal
			* sep, doublereal * work, integer * lwork, integer * iwork, integer * liwork, integer * info, ftnlen job_len, ftnlen compq_len)
{
	/* System generated locals */
	integer q_dim1, q_offset, t_dim1, t_offset, i__1, i__2;
	doublereal d__1, d__2;

	/* Builtin functions */
	double sqrt (doublereal);

	/* Local variables */
	static integer k, n1, n2, kk, nn, ks;
	static doublereal est;
	static integer kase;
	static logical pair;
	static integer ierr;
	static logical swap;
	static doublereal scale;
	static integer lwmin;
	static logical wantq, wants;
	static doublereal rnorm;
	static logical wantbh;
	static integer liwmin;
	static logical wantsp, lquery;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     June 30, 1999 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DTRSEN reorders the real Schur factorization of a real matrix */
	/*  A = Q*T*Q**T, so that a selected cluster of eigenvalues appears in */
	/*  the leading diagonal blocks of the upper quasi-triangular matrix T, */
	/*  and the leading columns of Q form an orthonormal basis of the */
	/*  corresponding right invariant subspace. */

	/*  Optionally the routine computes the reciprocal condition numbers of */
	/*  the cluster of eigenvalues and/or the invariant subspace. */

	/*  T must be in Schur canonical form (as returned by DHSEQR), that is, */
	/*  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each */
	/*  2-by-2 diagonal block has its diagonal elemnts equal and its */
	/*  off-diagonal elements of opposite sign. */

	/*  Arguments */
	/*  ========= */

	/*  JOB     (input) CHARACTER*1 */
	/*          Specifies whether condition numbers are required for the */
	/*          cluster of eigenvalues (S) or the invariant subspace (SEP): */
	/*          = 'N': none; */
	/*          = 'E': for eigenvalues only (S); */
	/*          = 'V': for invariant subspace only (SEP); */
	/*          = 'B': for both eigenvalues and invariant subspace (S and */
	/*                 SEP). */

	/*  COMPQ   (input) CHARACTER*1 */
	/*          = 'V': update the matrix Q of Schur vectors; */
	/*          = 'N': do not update Q. */

	/*  SELECT  (input) LOGICAL array, dimension (N) */
	/*          SELECT specifies the eigenvalues in the selected cluster. To */
	/*          select a real eigenvalue w(j), SELECT(j) must be set to */
	/*          .TRUE.. To select a complex conjugate pair of eigenvalues */
	/*          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block, */
	/*          either SELECT(j) or SELECT(j+1) or both must be set to */
	/*          .TRUE.; a complex conjugate pair of eigenvalues must be */
	/*          either both included in the cluster or both excluded. */

	/*  N       (input) INTEGER */
	/*          The order of the matrix T. N >= 0. */

	/*  T       (input/output) DOUBLE PRECISION array, dimension (LDT,N) */
	/*          On entry, the upper quasi-triangular matrix T, in Schur */
	/*          canonical form. */
	/*          On exit, T is overwritten by the reordered matrix T, again in */
	/*          Schur canonical form, with the selected eigenvalues in the */
	/*          leading diagonal blocks. */

	/*  LDT     (input) INTEGER */
	/*          The leading dimension of the array T. LDT >= max(1,N). */

	/*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
	/*          On entry, if COMPQ = 'V', the matrix Q of Schur vectors. */
	/*          On exit, if COMPQ = 'V', Q has been postmultiplied by the */
	/*          orthogonal transformation matrix which reorders T; the */
	/*          leading M columns of Q form an orthonormal basis for the */
	/*          specified invariant subspace. */
	/*          If COMPQ = 'N', Q is not referenced. */

	/*  LDQ     (input) INTEGER */
	/*          The leading dimension of the array Q. */
	/*          LDQ >= 1; and if COMPQ = 'V', LDQ >= N. */

	/*  WR      (output) DOUBLE PRECISION array, dimension (N) */
	/*  WI      (output) DOUBLE PRECISION array, dimension (N) */
	/*          The real and imaginary parts, respectively, of the reordered */
	/*          eigenvalues of T. The eigenvalues are stored in the same */
	/*          order as on the diagonal of T, with WR(i) = T(i,i) and, if */
	/*          T(i:i+1,i:i+1) is a 2-by-2 diagonal block, WI(i) > 0 and */
	/*          WI(i+1) = -WI(i). Note that if a complex eigenvalue is */
	/*          sufficiently ill-conditioned, then its value may differ */
	/*          significantly from its value before reordering. */

	/*  M       (output) INTEGER */
	/*          The dimension of the specified invariant subspace. */
	/*          0 < = M <= N. */

	/*  S       (output) DOUBLE PRECISION */
	/*          If JOB = 'E' or 'B', S is a lower bound on the reciprocal */
	/*          condition number for the selected cluster of eigenvalues. */
	/*          S cannot underestimate the true reciprocal condition number */
	/*          by more than a factor of sqrt(N). If M = 0 or N, S = 1. */
	/*          If JOB = 'N' or 'V', S is not referenced. */

	/*  SEP     (output) DOUBLE PRECISION */
	/*          If JOB = 'V' or 'B', SEP is the estimated reciprocal */
	/*          condition number of the specified invariant subspace. If */
	/*          M = 0 or N, SEP = norm(T). */
	/*          If JOB = 'N' or 'E', SEP is not referenced. */

	/*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK) */
	/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

	/*  LWORK   (input) INTEGER */
	/*          The dimension of the array WORK. */
	/*          If JOB = 'N', LWORK >= max(1,N); */
	/*          if JOB = 'E', LWORK >= M*(N-M); */
	/*          if JOB = 'V' or 'B', LWORK >= 2*M*(N-M). */

	/*          If LWORK = -1, then a workspace query is assumed; the routine */
	/*          only calculates the optimal size of the WORK array, returns */
	/*          this value as the first entry of the WORK array, and no error */
	/*          message related to LWORK is issued by XERBLA. */

	/*  IWORK   (workspace) INTEGER array, dimension (LIWORK) */
	/*          IF JOB = 'N' or 'E', IWORK is not referenced. */

	/*  LIWORK  (input) INTEGER */
	/*          The dimension of the array IWORK. */
	/*          If JOB = 'N' or 'E', LIWORK >= 1; */
	/*          if JOB = 'V' or 'B', LIWORK >= M*(N-M). */

	/*          If LIWORK = -1, then a workspace query is assumed; the */
	/*          routine only calculates the optimal size of the IWORK array, */
	/*          returns this value as the first entry of the IWORK array, and */
	/*          no error message related to LIWORK is issued by XERBLA. */

	/*  INFO    (output) INTEGER */
	/*          = 0: successful exit */
	/*          < 0: if INFO = -i, the i-th argument had an illegal value */
	/*          = 1: reordering of T failed because some eigenvalues are too */
	/*               close to separate (the problem is very ill-conditioned); */
	/*               T may have been partially reordered, and WR and WI */
	/*               contain the eigenvalues in the same order as in T; S and */
	/*               SEP (if requested) are set to zero. */

	/*  Further Details */
	/*  =============== */

	/*  DTRSEN first collects the selected eigenvalues by computing an */
	/*  orthogonal transformation Z to move them to the top left corner of T. */
	/*  In other words, the selected eigenvalues are the eigenvalues of T11 */
	/*  in: */

	/*                Z'*T*Z = ( T11 T12 ) n1 */
	/*                         (  0  T22 ) n2 */
	/*                            n1  n2 */

	/*  where N = n1+n2 and Z' means the transpose of Z. The first n1 columns */
	/*  of Z span the specified invariant subspace of T. */

	/*  If T has been obtained from the real Schur factorization of a matrix */
	/*  A = Q*T*Q', then the reordered real Schur factorization of A is given */
	/*  by A = (Q*Z)*(Z'*T*Z)*(Q*Z)', and the first n1 columns of Q*Z span */
	/*  the corresponding invariant subspace of A. */

	/*  The reciprocal condition number of the average of the eigenvalues of */
	/*  T11 may be returned in S. S lies between 0 (very badly conditioned) */
	/*  and 1 (very well conditioned). It is computed as follows. First we */
	/*  compute R so that */

	/*                         P = ( I  R ) n1 */
	/*                             ( 0  0 ) n2 */
	/*                               n1 n2 */

	/*  is the projector on the invariant subspace associated with T11. */
	/*  R is the solution of the Sylvester equation: */

	/*                        T11*R - R*T22 = T12. */

	/*  Let F-norm(M) denote the Frobenius-norm of M and 2-norm(M) denote */
	/*  the two-norm of M. Then S is computed as the lower bound */

	/*                      (1 + F-norm(R)**2)**(-1/2) */

	/*  on the reciprocal of 2-norm(P), the true reciprocal condition number. */
	/*  S cannot underestimate 1 / 2-norm(P) by more than a factor of */
	/*  sqrt(N). */

	/*  An approximate error bound for the computed average of the */
	/*  eigenvalues of T11 is */

	/*                         EPS * norm(T) / S */

	/*  where EPS is the machine precision. */

	/*  The reciprocal condition number of the right invariant subspace */
	/*  spanned by the first n1 columns of Z (or of Q*Z) is returned in SEP. */
	/*  SEP is defined as the separation of T11 and T22: */

	/*                     sep( T11, T22 ) = sigma-min( C ) */

	/*  where sigma-min(C) is the smallest singular value of the */
	/*  n1*n2-by-n1*n2 matrix */

	/*     C  = kprod( I(n2), T11 ) - kprod( transpose(T22), I(n1) ) */

	/*  I(m) is an m by m identity matrix, and kprod denotes the Kronecker */
	/*  product. We estimate sigma-min(C) by the reciprocal of an estimate of */
	/*  the 1-norm of inverse(C). The true reciprocal 1-norm of inverse(C) */
	/*  cannot differ from sigma-min(C) by more than a factor of sqrt(n1*n2). */

	/*  When SEP is small, small changes in T can cause large changes in */
	/*  the invariant subspace. An approximate bound on the maximum angular */
	/*  error in the computed right invariant subspace is */

	/*                      EPS * norm(T) / SEP */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Decode and test the input parameters */

	/* Parameter adjustments */
	--select;
	t_dim1 = *ldt;
	t_offset = 1 + t_dim1;
	t -= t_offset;
	q_dim1 = *ldq;
	q_offset = 1 + q_dim1;
	q -= q_offset;
	--wr;
	--wi;
	--work;
	--iwork;

	/* Function Body */
	wantbh = lsame_ (job, "B", (ftnlen) 1, (ftnlen) 1);
	wants = lsame_ (job, "E", (ftnlen) 1, (ftnlen) 1) || wantbh;
	wantsp = lsame_ (job, "V", (ftnlen) 1, (ftnlen) 1) || wantbh;
	wantq = lsame_ (compq, "V", (ftnlen) 1, (ftnlen) 1);

	*info = 0;
	lquery = *lwork == -1;
	if (!lsame_ (job, "N", (ftnlen) 1, (ftnlen) 1) && !wants && !wantsp)
	{
		*info = -1;
	}
	else if (!lsame_ (compq, "N", (ftnlen) 1, (ftnlen) 1) && !wantq)
	{
		*info = -2;
	}
	else if (*n < 0)
	{
		*info = -4;
	}
	else if (*ldt < max (1, *n))
	{
		*info = -6;
	}
	else if (*ldq < 1 || wantq && *ldq < *n)
	{
		*info = -8;
	}
	else
	{

		/*        Set M to the dimension of the specified invariant subspace, */
		/*        and test LWORK and LIWORK. */

		*m = 0;
		pair = FALSE_;
		i__1 = *n;
		for (k = 1; k <= i__1; ++k)
		{
			if (pair)
			{
				pair = FALSE_;
			}
			else
			{
				if (k < *n)
				{
					if (t[k + 1 + k * t_dim1] == 0.)
					{
						if (select[k])
						{
							++(*m);
						}
					}
					else
					{
						pair = TRUE_;
						if (select[k] || select[k + 1])
						{
							*m += 2;
						}
					}
				}
				else
				{
					if (select[*n])
					{
						++(*m);
					}
				}
			}
			/* L10: */
		}

		n1 = *m;
		n2 = *n - *m;
		nn = n1 * n2;

		if (wantsp)
		{
			/* Computing MAX */
			i__1 = 1, i__2 = nn << 1;
			lwmin = max (i__1, i__2);
			liwmin = max (1, nn);
		}
		else if (lsame_ (job, "N", (ftnlen) 1, (ftnlen) 1))
		{
			lwmin = max (1, *n);
			liwmin = 1;
		}
		else if (lsame_ (job, "E", (ftnlen) 1, (ftnlen) 1))
		{
			lwmin = max (1, nn);
			liwmin = 1;
		}

		if (*lwork < lwmin && !lquery)
		{
			*info = -15;
		}
		else if (*liwork < liwmin && !lquery)
		{
			*info = -17;
		}
	}

	if (*info == 0)
	{
		work[1] = (doublereal) lwmin;
		iwork[1] = liwmin;
	}

	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DTRSEN", &i__1, (ftnlen) 6);
		return 0;
	}
	else if (lquery)
	{
		return 0;
	}

	/*     Quick return if possible. */

	if (*m == *n || *m == 0)
	{
		if (wants)
		{
			*s = 1.;
		}
		if (wantsp)
		{
			*sep = dlange_ ("1", n, n, &t[t_offset], ldt, &work[1], (ftnlen) 1);
		}
		goto L40;
	}

	/*     Collect the selected blocks at the top-left corner of T. */

	ks = 0;
	pair = FALSE_;
	i__1 = *n;
	for (k = 1; k <= i__1; ++k)
	{
		if (pair)
		{
			pair = FALSE_;
		}
		else
		{
			swap = select[k];
			if (k < *n)
			{
				if (t[k + 1 + k * t_dim1] != 0.)
				{
					pair = TRUE_;
					swap = swap || select[k + 1];
				}
			}
			if (swap)
			{
				++ks;

				/*              Swap the K-th block to position KS. */

				ierr = 0;
				kk = k;
				if (k != ks)
				{
					dtrexc_ (compq, n, &t[t_offset], ldt, &q[q_offset], ldq, &kk, &ks, &work[1], &ierr, (ftnlen) 1);
				}
				if (ierr == 1 || ierr == 2)
				{

					/*                 Blocks too close to swap: exit. */

					*info = 1;
					if (wants)
					{
						*s = 0.;
					}
					if (wantsp)
					{
						*sep = 0.;
					}
					goto L40;
				}
				if (pair)
				{
					++ks;
				}
			}
		}
		/* L20: */
	}

	if (wants)
	{

		/*        Solve Sylvester equation for R: */

		/*           T11*R - R*T22 = scale*T12 */

		dlacpy_ ("F", &n1, &n2, &t[(n1 + 1) * t_dim1 + 1], ldt, &work[1], &n1, (ftnlen) 1);
		dtrsyl_ ("N", "N", &c_n1, &n1, &n2, &t[t_offset], ldt, &t[n1 + 1 + (n1 + 1) * t_dim1], ldt, &work[1], &n1, &scale, &ierr, (ftnlen) 1, (ftnlen) 1);

		/*        Estimate the reciprocal of the condition number of the cluster */
		/*        of eigenvalues. */

		rnorm = dlange_ ("F", &n1, &n2, &work[1], &n1, &work[1], (ftnlen) 1);
		if (rnorm == 0.)
		{
			*s = 1.;
		}
		else
		{
			*s = scale / (sqrt (scale * scale / rnorm + rnorm) * sqrt (rnorm));
		}
	}

	if (wantsp)
	{

		/*        Estimate sep(T11,T22). */

		est = 0.;
		kase = 0;
	 L30:
		dlacon_ (&nn, &work[nn + 1], &work[1], &iwork[1], &est, &kase);
		if (kase != 0)
		{
			if (kase == 1)
			{

				/*              Solve  T11*R - R*T22 = scale*X. */

				dtrsyl_ ("N", "N", &c_n1, &n1, &n2, &t[t_offset], ldt, &t[n1 + 1 + (n1 + 1) * t_dim1], ldt, &work[1], &n1, &scale, &ierr, (ftnlen) 1, (ftnlen) 1);
			}
			else
			{

				/*              Solve  T11'*R - R*T22' = scale*X. */

				dtrsyl_ ("T", "T", &c_n1, &n1, &n2, &t[t_offset], ldt, &t[n1 + 1 + (n1 + 1) * t_dim1], ldt, &work[1], &n1, &scale, &ierr, (ftnlen) 1, (ftnlen) 1);
			}
			goto L30;
		}

		*sep = scale / est;
	}

 L40:

	/*     Store the output eigenvalues in WR and WI. */

	i__1 = *n;
	for (k = 1; k <= i__1; ++k)
	{
		wr[k] = t[k + k * t_dim1];
		wi[k] = 0.;
		/* L50: */
	}
	i__1 = *n - 1;
	for (k = 1; k <= i__1; ++k)
	{
		if (t[k + 1 + k * t_dim1] != 0.)
		{
			wi[k] = sqrt ((d__1 = t[k + (k + 1) * t_dim1], abs (d__1))) * sqrt ((d__2 = t[k + 1 + k * t_dim1], abs (d__2)));
			wi[k + 1] = -wi[k];
		}
		/* L60: */
	}

	work[1] = (doublereal) lwmin;
	iwork[1] = liwmin;

	return 0;

	/*     End of DTRSEN */

}	/* dtrsen_ */

static int
dtrsna_ (char *job, char *howmny, logical * select,
			integer * n, doublereal * t, integer * ldt, doublereal * vl, integer *
			ldvl, doublereal * vr, integer * ldvr, doublereal * s, doublereal * sep,
			integer * mm, integer * m, doublereal * work, integer * ldwork, integer * iwork, integer * info, ftnlen job_len, ftnlen howmny_len)
{
	/* System generated locals */
	integer t_dim1, t_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, work_dim1, work_offset, i__1, i__2;
	doublereal d__1, d__2;

	/* Builtin functions */
	double sqrt (doublereal);

	/* Local variables */
	static integer i__, j, k, n2;
	static doublereal cs;
	static integer nn, ks;
	static doublereal sn, mu, eps, est;
	static integer kase;
	static doublereal cond;
	static logical pair;
	static integer ierr;
	static doublereal dumm, prod;
	static integer ifst;
	static doublereal lnrm;
	static integer ilst;
	static doublereal rnrm;
	static doublereal prod1, prod2, scale, delta;
	static logical wants;
	static doublereal dummy[1];
	static doublereal bignum;
	static logical wantbh;
	static logical somcon;
	static doublereal smlnum;
	static logical wantsp;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     September 30, 1994 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DTRSNA estimates reciprocal condition numbers for specified */
	/*  eigenvalues and/or right eigenvectors of a real upper */
	/*  quasi-triangular matrix T (or of any matrix Q*T*Q**T with Q */
	/*  orthogonal). */

	/*  T must be in Schur canonical form (as returned by DHSEQR), that is, */
	/*  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each */
	/*  2-by-2 diagonal block has its diagonal elements equal and its */
	/*  off-diagonal elements of opposite sign. */

	/*  Arguments */
	/*  ========= */

	/*  JOB     (input) CHARACTER*1 */
	/*          Specifies whether condition numbers are required for */
	/*          eigenvalues (S) or eigenvectors (SEP): */
	/*          = 'E': for eigenvalues only (S); */
	/*          = 'V': for eigenvectors only (SEP); */
	/*          = 'B': for both eigenvalues and eigenvectors (S and SEP). */

	/*  HOWMNY  (input) CHARACTER*1 */
	/*          = 'A': compute condition numbers for all eigenpairs; */
	/*          = 'S': compute condition numbers for selected eigenpairs */
	/*                 specified by the array SELECT. */

	/*  SELECT  (input) LOGICAL array, dimension (N) */
	/*          If HOWMNY = 'S', SELECT specifies the eigenpairs for which */
	/*          condition numbers are required. To select condition numbers */
	/*          for the eigenpair corresponding to a real eigenvalue w(j), */
	/*          SELECT(j) must be set to .TRUE.. To select condition numbers */
	/*          corresponding to a complex conjugate pair of eigenvalues w(j) */
	/*          and w(j+1), either SELECT(j) or SELECT(j+1) or both, must be */
	/*          set to .TRUE.. */
	/*          If HOWMNY = 'A', SELECT is not referenced. */

	/*  N       (input) INTEGER */
	/*          The order of the matrix T. N >= 0. */

	/*  T       (input) DOUBLE PRECISION array, dimension (LDT,N) */
	/*          The upper quasi-triangular matrix T, in Schur canonical form. */

	/*  LDT     (input) INTEGER */
	/*          The leading dimension of the array T. LDT >= max(1,N). */

	/*  VL      (input) DOUBLE PRECISION array, dimension (LDVL,M) */
	/*          If JOB = 'E' or 'B', VL must contain left eigenvectors of T */
	/*          (or of any Q*T*Q**T with Q orthogonal), corresponding to the */
	/*          eigenpairs specified by HOWMNY and SELECT. The eigenvectors */
	/*          must be stored in consecutive columns of VL, as returned by */
	/*          DHSEIN or DTREVC. */
	/*          If JOB = 'V', VL is not referenced. */

	/*  LDVL    (input) INTEGER */
	/*          The leading dimension of the array VL. */
	/*          LDVL >= 1; and if JOB = 'E' or 'B', LDVL >= N. */

	/*  VR      (input) DOUBLE PRECISION array, dimension (LDVR,M) */
	/*          If JOB = 'E' or 'B', VR must contain right eigenvectors of T */
	/*          (or of any Q*T*Q**T with Q orthogonal), corresponding to the */
	/*          eigenpairs specified by HOWMNY and SELECT. The eigenvectors */
	/*          must be stored in consecutive columns of VR, as returned by */
	/*          DHSEIN or DTREVC. */
	/*          If JOB = 'V', VR is not referenced. */

	/*  LDVR    (input) INTEGER */
	/*          The leading dimension of the array VR. */
	/*          LDVR >= 1; and if JOB = 'E' or 'B', LDVR >= N. */

	/*  S       (output) DOUBLE PRECISION array, dimension (MM) */
	/*          If JOB = 'E' or 'B', the reciprocal condition numbers of the */
	/*          selected eigenvalues, stored in consecutive elements of the */
	/*          array. For a complex conjugate pair of eigenvalues two */
	/*          consecutive elements of S are set to the same value. Thus */
	/*          S(j), SEP(j), and the j-th columns of VL and VR all */
	/*          correspond to the same eigenpair (but not in general the */
	/*          j-th eigenpair, unless all eigenpairs are selected). */
	/*          If JOB = 'V', S is not referenced. */

	/*  SEP     (output) DOUBLE PRECISION array, dimension (MM) */
	/*          If JOB = 'V' or 'B', the estimated reciprocal condition */
	/*          numbers of the selected eigenvectors, stored in consecutive */
	/*          elements of the array. For a complex eigenvector two */
	/*          consecutive elements of SEP are set to the same value. If */
	/*          the eigenvalues cannot be reordered to compute SEP(j), SEP(j) */
	/*          is set to 0; this can only occur when the true value would be */
	/*          very small anyway. */
	/*          If JOB = 'E', SEP is not referenced. */

	/*  MM      (input) INTEGER */
	/*          The number of elements in the arrays S (if JOB = 'E' or 'B') */
	/*           and/or SEP (if JOB = 'V' or 'B'). MM >= M. */

	/*  M       (output) INTEGER */
	/*          The number of elements of the arrays S and/or SEP actually */
	/*          used to store the estimated condition numbers. */
	/*          If HOWMNY = 'A', M is set to N. */

	/*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,N+1) */
	/*          If JOB = 'E', WORK is not referenced. */

	/*  LDWORK  (input) INTEGER */
	/*          The leading dimension of the array WORK. */
	/*          LDWORK >= 1; and if JOB = 'V' or 'B', LDWORK >= N. */

	/*  IWORK   (workspace) INTEGER array, dimension (N) */
	/*          If JOB = 'E', IWORK is not referenced. */

	/*  INFO    (output) INTEGER */
	/*          = 0: successful exit */
	/*          < 0: if INFO = -i, the i-th argument had an illegal value */

	/*  Further Details */
	/*  =============== */

	/*  The reciprocal of the condition number of an eigenvalue lambda is */
	/*  defined as */

	/*          S(lambda) = |v'*u| / (norm(u)*norm(v)) */

	/*  where u and v are the right and left eigenvectors of T corresponding */
	/*  to lambda; v' denotes the conjugate-transpose of v, and norm(u) */
	/*  denotes the Euclidean norm. These reciprocal condition numbers always */
	/*  lie between zero (very badly conditioned) and one (very well */
	/*  conditioned). If n = 1, S(lambda) is defined to be 1. */

	/*  An approximate error bound for a computed eigenvalue W(i) is given by */

	/*                      EPS * norm(T) / S(i) */

	/*  where EPS is the machine precision. */

	/*  The reciprocal of the condition number of the right eigenvector u */
	/*  corresponding to lambda is defined as follows. Suppose */

	/*              T = ( lambda  c  ) */
	/*                  (   0    T22 ) */

	/*  Then the reciprocal condition number is */

	/*          SEP( lambda, T22 ) = sigma-min( T22 - lambda*I ) */

	/*  where sigma-min denotes the smallest singular value. We approximate */
	/*  the smallest singular value by the reciprocal of an estimate of the */
	/*  one-norm of the inverse of T22 - lambda*I. If n = 1, SEP(1) is */
	/*  defined to be abs(T(1,1)). */

	/*  An approximate error bound for a computed right eigenvector VR(i) */
	/*  is given by */

	/*                      EPS * norm(T) / SEP(i) */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Local Arrays .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Decode and test the input parameters */

	/* Parameter adjustments */
	--select;
	t_dim1 = *ldt;
	t_offset = 1 + t_dim1;
	t -= t_offset;
	vl_dim1 = *ldvl;
	vl_offset = 1 + vl_dim1;
	vl -= vl_offset;
	vr_dim1 = *ldvr;
	vr_offset = 1 + vr_dim1;
	vr -= vr_offset;
	--s;
	--sep;
	work_dim1 = *ldwork;
	work_offset = 1 + work_dim1;
	work -= work_offset;
	--iwork;

	/* Function Body */
	wantbh = lsame_ (job, "B", (ftnlen) 1, (ftnlen) 1);
	wants = lsame_ (job, "E", (ftnlen) 1, (ftnlen) 1) || wantbh;
	wantsp = lsame_ (job, "V", (ftnlen) 1, (ftnlen) 1) || wantbh;

	somcon = lsame_ (howmny, "S", (ftnlen) 1, (ftnlen) 1);

	*info = 0;
	if (!wants && !wantsp)
	{
		*info = -1;
	}
	else if (!lsame_ (howmny, "A", (ftnlen) 1, (ftnlen) 1) && !somcon)
	{
		*info = -2;
	}
	else if (*n < 0)
	{
		*info = -4;
	}
	else if (*ldt < max (1, *n))
	{
		*info = -6;
	}
	else if (*ldvl < 1 || wants && *ldvl < *n)
	{
		*info = -8;
	}
	else if (*ldvr < 1 || wants && *ldvr < *n)
	{
		*info = -10;
	}
	else
	{

		/*        Set M to the number of eigenpairs for which condition numbers */
		/*        are required, and test MM. */

		if (somcon)
		{
			*m = 0;
			pair = FALSE_;
			i__1 = *n;
			for (k = 1; k <= i__1; ++k)
			{
				if (pair)
				{
					pair = FALSE_;
				}
				else
				{
					if (k < *n)
					{
						if (t[k + 1 + k * t_dim1] == 0.)
						{
							if (select[k])
							{
								++(*m);
							}
						}
						else
						{
							pair = TRUE_;
							if (select[k] || select[k + 1])
							{
								*m += 2;
							}
						}
					}
					else
					{
						if (select[*n])
						{
							++(*m);
						}
					}
				}
				/* L10: */
			}
		}
		else
		{
			*m = *n;
		}

		if (*mm < *m)
		{
			*info = -13;
		}
		else if (*ldwork < 1 || wantsp && *ldwork < *n)
		{
			*info = -16;
		}
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DTRSNA", &i__1, (ftnlen) 6);
		return 0;
	}

	/*     Quick return if possible */

	if (*n == 0)
	{
		return 0;
	}

	if (*n == 1)
	{
		if (somcon)
		{
			if (!select[1])
			{
				return 0;
			}
		}
		if (wants)
		{
			s[1] = 1.;
		}
		if (wantsp)
		{
			sep[1] = (d__1 = t[t_dim1 + 1], abs (d__1));
		}
		return 0;
	}

	/*     Get machine constants */

	eps = dlamch_ ("P", (ftnlen) 1);
	smlnum = dlamch_ ("S", (ftnlen) 1) / eps;
	bignum = 1. / smlnum;
	dlabad_ (&smlnum, &bignum);

	ks = 0;
	pair = FALSE_;
	i__1 = *n;
	for (k = 1; k <= i__1; ++k)
	{

		/*        Determine whether T(k,k) begins a 1-by-1 or 2-by-2 block. */

		if (pair)
		{
			pair = FALSE_;
			goto L60;
		}
		else
		{
			if (k < *n)
			{
				pair = t[k + 1 + k * t_dim1] != 0.;
			}
		}

		/*        Determine whether condition numbers are required for the k-th */
		/*        eigenpair. */

		if (somcon)
		{
			if (pair)
			{
				if (!select[k] && !select[k + 1])
				{
					goto L60;
				}
			}
			else
			{
				if (!select[k])
				{
					goto L60;
				}
			}
		}

		++ks;

		if (wants)
		{

			/*           Compute the reciprocal condition number of the k-th */
			/*           eigenvalue. */

			if (!pair)
			{

				/*              Real eigenvalue. */

				prod = ddot_ (n, &vr[ks * vr_dim1 + 1], &c__1, &vl[ks * vl_dim1 + 1], &c__1);
				rnrm = dnrm2_ (n, &vr[ks * vr_dim1 + 1], &c__1);
				lnrm = dnrm2_ (n, &vl[ks * vl_dim1 + 1], &c__1);
				s[ks] = abs (prod) / (rnrm * lnrm);
			}
			else
			{

				/*              Complex eigenvalue. */

				prod1 = ddot_ (n, &vr[ks * vr_dim1 + 1], &c__1, &vl[ks * vl_dim1 + 1], &c__1);
				prod1 += ddot_ (n, &vr[(ks + 1) * vr_dim1 + 1], &c__1, &vl[(ks + 1) * vl_dim1 + 1], &c__1);
				prod2 = ddot_ (n, &vl[ks * vl_dim1 + 1], &c__1, &vr[(ks + 1) * vr_dim1 + 1], &c__1);
				prod2 -= ddot_ (n, &vl[(ks + 1) * vl_dim1 + 1], &c__1, &vr[ks * vr_dim1 + 1], &c__1);
				d__1 = dnrm2_ (n, &vr[ks * vr_dim1 + 1], &c__1);
				d__2 = dnrm2_ (n, &vr[(ks + 1) * vr_dim1 + 1], &c__1);
				rnrm = dlapy2_ (&d__1, &d__2);
				d__1 = dnrm2_ (n, &vl[ks * vl_dim1 + 1], &c__1);
				d__2 = dnrm2_ (n, &vl[(ks + 1) * vl_dim1 + 1], &c__1);
				lnrm = dlapy2_ (&d__1, &d__2);
				cond = dlapy2_ (&prod1, &prod2) / (rnrm * lnrm);
				s[ks] = cond;
				s[ks + 1] = cond;
			}
		}

		if (wantsp)
		{

			/*           Estimate the reciprocal condition number of the k-th */
			/*           eigenvector. */

			/*           Copy the matrix T to the array WORK and swap the diagonal */
			/*           block beginning at T(k,k) to the (1,1) position. */

			dlacpy_ ("Full", n, n, &t[t_offset], ldt, &work[work_offset], ldwork, (ftnlen) 4);
			ifst = k;
			ilst = 1;
			dtrexc_ ("No Q", n, &work[work_offset], ldwork, dummy, &c__1, &ifst, &ilst, &work[(*n + 1) * work_dim1 + 1], &ierr, (ftnlen) 4);

			if (ierr == 1 || ierr == 2)
			{

				/*              Could not swap because blocks not well separated */

				scale = 1.;
				est = bignum;
			}
			else
			{

				/*              Reordering successful */

				if (work[work_dim1 + 2] == 0.)
				{

					/*                 Form C = T22 - lambda*I in WORK(2:N,2:N). */

					i__2 = *n;
					for (i__ = 2; i__ <= i__2; ++i__)
					{
						work[i__ + i__ * work_dim1] -= work[work_dim1 + 1];
						/* L20: */
					}
					n2 = 1;
					nn = *n - 1;
				}
				else
				{

					/*                 Triangularize the 2 by 2 block by unitary */
					/*                 transformation U = [  cs   i*ss ] */
					/*                                    [ i*ss   cs  ]. */
					/*                 such that the (1,1) position of WORK is complex */
					/*                 eigenvalue lambda with positive imaginary part. (2,2) */
					/*                 position of WORK is the complex eigenvalue lambda */
					/*                 with negative imaginary  part. */

					mu = sqrt ((d__1 = work[(work_dim1 << 1) + 1], abs (d__1))) * sqrt ((d__2 = work[work_dim1 + 2], abs (d__2)));
					delta = dlapy2_ (&mu, &work[work_dim1 + 2]);
					cs = mu / delta;
					sn = -work[work_dim1 + 2] / delta;

					/*                 Form */

					/*                 C' = WORK(2:N,2:N) + i*[rwork(1) ..... rwork(n-1) ] */
					/*                                        [   mu                     ] */
					/*                                        [         ..               ] */
					/*                                        [             ..           ] */
					/*                                        [                  mu      ] */
					/*                 where C' is conjugate transpose of complex matrix C, */
					/*                 and RWORK is stored starting in the N+1-st column of */
					/*                 WORK. */

					i__2 = *n;
					for (j = 3; j <= i__2; ++j)
					{
						work[j * work_dim1 + 2] = cs * work[j * work_dim1 + 2];
						work[j + j * work_dim1] -= work[work_dim1 + 1];
						/* L30: */
					}
					work[(work_dim1 << 1) + 2] = 0.;

					work[(*n + 1) * work_dim1 + 1] = mu * 2.;
					i__2 = *n - 1;
					for (i__ = 2; i__ <= i__2; ++i__)
					{
						work[i__ + (*n + 1) * work_dim1] = sn * work[(i__ + 1) * work_dim1 + 1];
						/* L40: */
					}
					n2 = 2;
					nn = *n - 1 << 1;
				}

				/*              Estimate norm(inv(C')) */

				est = 0.;
				kase = 0;
			 L50:
				dlacon_ (&nn, &work[(*n + 2) * work_dim1 + 1], &work[(*n + 4) * work_dim1 + 1], &iwork[1], &est, &kase);
				if (kase != 0)
				{
					if (kase == 1)
					{
						if (n2 == 1)
						{

							/*                       Real eigenvalue: solve C'*x = scale*c. */

							i__2 = *n - 1;
							dlaqtr_ (&c_true, &c_true, &i__2, &work[(work_dim1
																				  << 1) + 2], ldwork, dummy, &dumm, &scale,
										&work[(*n + 4) * work_dim1 + 1], &work[(*n + 6) * work_dim1 + 1], &ierr);
						}
						else
						{

							/*                       Complex eigenvalue: solve */
							/*                       C'*(p+iq) = scale*(c+id) in real arithmetic. */

							i__2 = *n - 1;
							dlaqtr_ (&c_true, &c_false, &i__2, &work[(work_dim1 << 1) + 2], ldwork, &work[(*n +
																																	 1) * work_dim1 + 1], &mu, &scale,
										&work[(*n + 4) * work_dim1 + 1], &work[(*n + 6) * work_dim1 + 1], &ierr);
						}
					}
					else
					{
						if (n2 == 1)
						{

							/*                       Real eigenvalue: solve C*x = scale*c. */

							i__2 = *n - 1;
							dlaqtr_ (&c_false, &c_true, &i__2, &work[(work_dim1 << 1) + 2], ldwork, dummy, &dumm, &scale, &work[(*n + 4) * work_dim1
																																								 + 1], &work[(*n + 6) * work_dim1 + 1],
										&ierr);
						}
						else
						{

							/*                       Complex eigenvalue: solve */
							/*                       C*(p+iq) = scale*(c+id) in real arithmetic. */

							i__2 = *n - 1;
							dlaqtr_ (&c_false, &c_false, &i__2, &work[(work_dim1 << 1) + 2], ldwork, &work[(*n +
																																	  1) * work_dim1 + 1], &mu, &scale,
										&work[(*n + 4) * work_dim1 + 1], &work[(*n + 6) * work_dim1 + 1], &ierr);

						}
					}

					goto L50;
				}
			}

			sep[ks] = scale / max (est, smlnum);
			if (pair)
			{
				sep[ks + 1] = sep[ks];
			}
		}

		if (pair)
		{
			++ks;
		}

	 L60:
		;
	}
	return 0;

	/*     End of DTRSNA */

}	/* dtrsna_ */

static int
dtrsyl_ (char *trana, char *tranb, integer * isgn, integer
			* m, integer * n, doublereal * a, integer * lda, doublereal * b, integer *
			ldb, doublereal * c__, integer * ldc, doublereal * scale, integer * info, ftnlen trana_len, ftnlen tranb_len)
{
	/* System generated locals */
	integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, i__3, i__4;
	doublereal d__1, d__2;

	/* Local variables */
	static integer j, k, l;
	static doublereal x[4] /* was [2][2] */ ;
	static integer k1, k2, l1, l2;
	static doublereal a11, db, da11, vec[4] /* was [2][2] */ , dum[1], eps, sgn;
	static integer ierr;
	static doublereal smin, suml, sumr;
	static integer knext, lnext;
	static doublereal xnorm;
	static doublereal scaloc;
	static doublereal bignum;
	static logical notrna, notrnb;
	static doublereal smlnum;


	/*  -- LAPACK routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     March 31, 1993 */

	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  DTRSYL solves the real Sylvester matrix equation: */

	/*     op(A)*X + X*op(B) = scale*C or */
	/*     op(A)*X - X*op(B) = scale*C, */

	/*  where op(A) = A or A**T, and  A and B are both upper quasi- */
	/*  triangular. A is M-by-M and B is N-by-N; the right hand side C and */
	/*  the solution X are M-by-N; and scale is an output scale factor, set */
	/*  <= 1 to avoid overflow in X. */

	/*  A and B must be in Schur canonical form (as returned by DHSEQR), that */
	/*  is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; */
	/*  each 2-by-2 diagonal block has its diagonal elements equal and its */
	/*  off-diagonal elements of opposite sign. */

	/*  Arguments */
	/*  ========= */

	/*  TRANA   (input) CHARACTER*1 */
	/*          Specifies the option op(A): */
	/*          = 'N': op(A) = A    (No transpose) */
	/*          = 'T': op(A) = A**T (Transpose) */
	/*          = 'C': op(A) = A**H (Conjugate transpose = Transpose) */

	/*  TRANB   (input) CHARACTER*1 */
	/*          Specifies the option op(B): */
	/*          = 'N': op(B) = B    (No transpose) */
	/*          = 'T': op(B) = B**T (Transpose) */
	/*          = 'C': op(B) = B**H (Conjugate transpose = Transpose) */

	/*  ISGN    (input) INTEGER */
	/*          Specifies the sign in the equation: */
	/*          = +1: solve op(A)*X + X*op(B) = scale*C */
	/*          = -1: solve op(A)*X - X*op(B) = scale*C */

	/*  M       (input) INTEGER */
	/*          The order of the matrix A, and the number of rows in the */
	/*          matrices X and C. M >= 0. */

	/*  N       (input) INTEGER */
	/*          The order of the matrix B, and the number of columns in the */
	/*          matrices X and C. N >= 0. */

	/*  A       (input) DOUBLE PRECISION array, dimension (LDA,M) */
	/*          The upper quasi-triangular matrix A, in Schur canonical form. */

	/*  LDA     (input) INTEGER */
	/*          The leading dimension of the array A. LDA >= max(1,M). */

	/*  B       (input) DOUBLE PRECISION array, dimension (LDB,N) */
	/*          The upper quasi-triangular matrix B, in Schur canonical form. */

	/*  LDB     (input) INTEGER */
	/*          The leading dimension of the array B. LDB >= max(1,N). */

	/*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
	/*          On entry, the M-by-N right hand side matrix C. */
	/*          On exit, C is overwritten by the solution matrix X. */

	/*  LDC     (input) INTEGER */
	/*          The leading dimension of the array C. LDC >= max(1,M) */

	/*  SCALE   (output) DOUBLE PRECISION */
	/*          The scale factor, scale, set <= 1 to avoid overflow in X. */

	/*  INFO    (output) INTEGER */
	/*          = 0: successful exit */
	/*          < 0: if INFO = -i, the i-th argument had an illegal value */
	/*          = 1: A and B have common or very close eigenvalues; perturbed */
	/*               values were used to solve the equation (but the matrices */
	/*               A and B are unchanged). */

	/*  ===================================================================== */

	/*     .. Parameters .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Local Arrays .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. External Subroutines .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	/*     Decode and Test input parameters */

	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;
	b_dim1 = *ldb;
	b_offset = 1 + b_dim1;
	b -= b_offset;
	c_dim1 = *ldc;
	c_offset = 1 + c_dim1;
	c__ -= c_offset;

	/* Function Body */
	notrna = lsame_ (trana, "N", (ftnlen) 1, (ftnlen) 1);
	notrnb = lsame_ (tranb, "N", (ftnlen) 1, (ftnlen) 1);

	*info = 0;
	if (!notrna && !lsame_ (trana, "T", (ftnlen) 1, (ftnlen) 1) && !lsame_ (trana, "C", (ftnlen) 1, (ftnlen) 1))
	{
		*info = -1;
	}
	else if (!notrnb && !lsame_ (tranb, "T", (ftnlen) 1, (ftnlen) 1) && !lsame_ (tranb, "C", (ftnlen) 1, (ftnlen) 1))
	{
		*info = -2;
	}
	else if (*isgn != 1 && *isgn != -1)
	{
		*info = -3;
	}
	else if (*m < 0)
	{
		*info = -4;
	}
	else if (*n < 0)
	{
		*info = -5;
	}
	else if (*lda < max (1, *m))
	{
		*info = -7;
	}
	else if (*ldb < max (1, *n))
	{
		*info = -9;
	}
	else if (*ldc < max (1, *m))
	{
		*info = -11;
	}
	if (*info != 0)
	{
		i__1 = -(*info);
		xerbla_ ("DTRSYL", &i__1, (ftnlen) 6);
		return 0;
	}

	/*     Quick return if possible */

	if (*m == 0 || *n == 0)
	{
		return 0;
	}

	/*     Set constants to control overflow */

	eps = dlamch_ ("P", (ftnlen) 1);
	smlnum = dlamch_ ("S", (ftnlen) 1);
	bignum = 1. / smlnum;
	dlabad_ (&smlnum, &bignum);
	smlnum = smlnum * (doublereal) (*m * *n) / eps;
	bignum = 1. / smlnum;

	/* Computing MAX */
	d__1 = smlnum, d__2 = eps * dlange_ ("M", m, m, &a[a_offset], lda, dum, (ftnlen) 1), d__1 = max (d__1, d__2), d__2 = eps * dlange_ ("M", n, n,
																																													&b[b_offset], ldb, dum,
																																													(ftnlen) 1);
	smin = max (d__1, d__2);

	*scale = 1.;
	sgn = (doublereal) (*isgn);

	if (notrna && notrnb)
	{

		/*        Solve    A*X + ISGN*X*B = scale*C. */

		/*        The (K,L)th block of X is determined starting from */
		/*        bottom-left corner column by column by */

		/*         A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L) */

		/*        Where */
		/*                  M                         L-1 */
		/*        R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)]. */
		/*                I=K+1                       J=1 */

		/*        Start column loop (index = L) */
		/*        L1 (L2) : column index of the first (first) row of X(K,L). */

		lnext = 1;
		i__1 = *n;
		for (l = 1; l <= i__1; ++l)
		{
			if (l < lnext)
			{
				goto L60;
			}
			if (l == *n)
			{
				l1 = l;
				l2 = l;
			}
			else
			{
				if (b[l + 1 + l * b_dim1] != 0.)
				{
					l1 = l;
					l2 = l + 1;
					lnext = l + 2;
				}
				else
				{
					l1 = l;
					l2 = l;
					lnext = l + 1;
				}
			}

			/*           Start row loop (index = K) */
			/*           K1 (K2): row index of the first (last) row of X(K,L). */

			knext = *m;
			for (k = *m; k >= 1; --k)
			{
				if (k > knext)
				{
					goto L50;
				}
				if (k == 1)
				{
					k1 = k;
					k2 = k;
				}
				else
				{
					if (a[k + (k - 1) * a_dim1] != 0.)
					{
						k1 = k - 1;
						k2 = k;
						knext = k - 2;
					}
					else
					{
						k1 = k;
						k2 = k;
						knext = k - 1;
					}
				}

				if (l1 == l2 && k1 == k2)
				{
					i__2 = *m - k1;
					/* Computing MIN */
					i__3 = k1 + 1;
					/* Computing MIN */
					i__4 = k1 + 1;
					suml = ddot_ (&i__2, &a[k1 + min (i__3, *m) * a_dim1], lda, &c__[min (i__4, *m) + l1 * c_dim1], &c__1);
					i__2 = l1 - 1;
					sumr = ddot_ (&i__2, &c__[k1 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
					vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
					scaloc = 1.;

					a11 = a[k1 + k1 * a_dim1] + sgn * b[l1 + l1 * b_dim1];
					da11 = abs (a11);
					if (da11 <= smin)
					{
						a11 = smin;
						da11 = smin;
						*info = 1;
					}
					db = abs (vec[0]);
					if (da11 < 1. && db > 1.)
					{
						if (db > bignum * da11)
						{
							scaloc = 1. / db;
						}
					}
					x[0] = vec[0] * scaloc / a11;

					if (scaloc != 1.)
					{
						i__2 = *n;
						for (j = 1; j <= i__2; ++j)
						{
							dscal_ (m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
							/* L10: */
						}
						*scale *= scaloc;
					}
					c__[k1 + l1 * c_dim1] = x[0];

				}
				else if (l1 == l2 && k1 != k2)
				{

					i__2 = *m - k2;
					/* Computing MIN */
					i__3 = k2 + 1;
					/* Computing MIN */
					i__4 = k2 + 1;
					suml = ddot_ (&i__2, &a[k1 + min (i__3, *m) * a_dim1], lda, &c__[min (i__4, *m) + l1 * c_dim1], &c__1);
					i__2 = l1 - 1;
					sumr = ddot_ (&i__2, &c__[k1 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
					vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

					i__2 = *m - k2;
					/* Computing MIN */
					i__3 = k2 + 1;
					/* Computing MIN */
					i__4 = k2 + 1;
					suml = ddot_ (&i__2, &a[k2 + min (i__3, *m) * a_dim1], lda, &c__[min (i__4, *m) + l1 * c_dim1], &c__1);
					i__2 = l1 - 1;
					sumr = ddot_ (&i__2, &c__[k2 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
					vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

					d__1 = -sgn * b[l1 + l1 * b_dim1];
					dlaln2_ (&c_false, &c__2, &c__1, &smin, &c_b348, &a[k1 +
																						 k1 * a_dim1], lda, &c_b348, &c_b348, vec, &c__2, &d__1, &c_b507, x, &c__2, &scaloc, &xnorm,
								&ierr);
					if (ierr != 0)
					{
						*info = 1;
					}

					if (scaloc != 1.)
					{
						i__2 = *n;
						for (j = 1; j <= i__2; ++j)
						{
							dscal_ (m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
							/* L20: */
						}
						*scale *= scaloc;
					}
					c__[k1 + l1 * c_dim1] = x[0];
					c__[k2 + l1 * c_dim1] = x[1];

				}
				else if (l1 != l2 && k1 == k2)
				{

					i__2 = *m - k1;
					/* Computing MIN */
					i__3 = k1 + 1;
					/* Computing MIN */
					i__4 = k1 + 1;
					suml = ddot_ (&i__2, &a[k1 + min (i__3, *m) * a_dim1], lda, &c__[min (i__4, *m) + l1 * c_dim1], &c__1);
					i__2 = l1 - 1;
					sumr = ddot_ (&i__2, &c__[k1 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
					vec[0] = sgn * (c__[k1 + l1 * c_dim1] - (suml + sgn * sumr));

					i__2 = *m - k1;
					/* Computing MIN */
					i__3 = k1 + 1;
					/* Computing MIN */
					i__4 = k1 + 1;
					suml = ddot_ (&i__2, &a[k1 + min (i__3, *m) * a_dim1], lda, &c__[min (i__4, *m) + l2 * c_dim1], &c__1);
					i__2 = l1 - 1;
					sumr = ddot_ (&i__2, &c__[k1 + c_dim1], ldc, &b[l2 * b_dim1 + 1], &c__1);
					vec[1] = sgn * (c__[k1 + l2 * c_dim1] - (suml + sgn * sumr));

					d__1 = -sgn * a[k1 + k1 * a_dim1];
					dlaln2_ (&c_true, &c__2, &c__1, &smin, &c_b348, &b[l1 + l1
																						* b_dim1], ldb, &c_b348, &c_b348, vec, &c__2, &d__1, &c_b507, x, &c__2, &scaloc, &xnorm,
								&ierr);
					if (ierr != 0)
					{
						*info = 1;
					}

					if (scaloc != 1.)
					{
						i__2 = *n;
						for (j = 1; j <= i__2; ++j)
						{
							dscal_ (m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
							/* L30: */
						}
						*scale *= scaloc;
					}
					c__[k1 + l1 * c_dim1] = x[0];
					c__[k1 + l2 * c_dim1] = x[1];

				}
				else if (l1 != l2 && k1 != k2)
				{

					i__2 = *m - k2;
					/* Computing MIN */
					i__3 = k2 + 1;
					/* Computing MIN */
					i__4 = k2 + 1;
					suml = ddot_ (&i__2, &a[k1 + min (i__3, *m) * a_dim1], lda, &c__[min (i__4, *m) + l1 * c_dim1], &c__1);
					i__2 = l1 - 1;
					sumr = ddot_ (&i__2, &c__[k1 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
					vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

					i__2 = *m - k2;
					/* Computing MIN */
					i__3 = k2 + 1;
					/* Computing MIN */
					i__4 = k2 + 1;
					suml = ddot_ (&i__2, &a[k1 + min (i__3, *m) * a_dim1], lda, &c__[min (i__4, *m) + l2 * c_dim1], &c__1);
					i__2 = l1 - 1;
					sumr = ddot_ (&i__2, &c__[k1 + c_dim1], ldc, &b[l2 * b_dim1 + 1], &c__1);
					vec[2] = c__[k1 + l2 * c_dim1] - (suml + sgn * sumr);

					i__2 = *m - k2;
					/* Computing MIN */
					i__3 = k2 + 1;
					/* Computing MIN */
					i__4 = k2 + 1;
					suml = ddot_ (&i__2, &a[k2 + min (i__3, *m) * a_dim1], lda, &c__[min (i__4, *m) + l1 * c_dim1], &c__1);
					i__2 = l1 - 1;
					sumr = ddot_ (&i__2, &c__[k2 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
					vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

					i__2 = *m - k2;
					/* Computing MIN */
					i__3 = k2 + 1;
					/* Computing MIN */
					i__4 = k2 + 1;
					suml = ddot_ (&i__2, &a[k2 + min (i__3, *m) * a_dim1], lda, &c__[min (i__4, *m) + l2 * c_dim1], &c__1);
					i__2 = l1 - 1;
					sumr = ddot_ (&i__2, &c__[k2 + c_dim1], ldc, &b[l2 * b_dim1 + 1], &c__1);
					vec[3] = c__[k2 + l2 * c_dim1] - (suml + sgn * sumr);

					dlasy2_ (&c_false, &c_false, isgn, &c__2, &c__2, &a[k1 +
																						 k1 * a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, vec, &c__2, &scaloc, x, &c__2, &xnorm, &ierr);
					if (ierr != 0)
					{
						*info = 1;
					}

					if (scaloc != 1.)
					{
						i__2 = *n;
						for (j = 1; j <= i__2; ++j)
						{
							dscal_ (m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
							/* L40: */
						}
						*scale *= scaloc;
					}
					c__[k1 + l1 * c_dim1] = x[0];
					c__[k1 + l2 * c_dim1] = x[2];
					c__[k2 + l1 * c_dim1] = x[1];
					c__[k2 + l2 * c_dim1] = x[3];
				}

			 L50:
				;
			}

		 L60:
			;
		}

	}
	else if (!notrna && notrnb)
	{

		/*        Solve    A' *X + ISGN*X*B = scale*C. */

		/*        The (K,L)th block of X is determined starting from */
		/*        upper-left corner column by column by */

		/*          A(K,K)'*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L) */

		/*        Where */
		/*                   K-1                        L-1 */
		/*          R(K,L) = SUM [A(I,K)'*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)] */
		/*                   I=1                        J=1 */

		/*        Start column loop (index = L) */
		/*        L1 (L2): column index of the first (last) row of X(K,L) */

		lnext = 1;
		i__1 = *n;
		for (l = 1; l <= i__1; ++l)
		{
			if (l < lnext)
			{
				goto L120;
			}
			if (l == *n)
			{
				l1 = l;
				l2 = l;
			}
			else
			{
				if (b[l + 1 + l * b_dim1] != 0.)
				{
					l1 = l;
					l2 = l + 1;
					lnext = l + 2;
				}
				else
				{
					l1 = l;
					l2 = l;
					lnext = l + 1;
				}
			}

			/*           Start row loop (index = K) */
			/*           K1 (K2): row index of the first (last) row of X(K,L) */

			knext = 1;
			i__2 = *m;
			for (k = 1; k <= i__2; ++k)
			{
				if (k < knext)
				{
					goto L110;
				}
				if (k == *m)
				{
					k1 = k;
					k2 = k;
				}
				else
				{
					if (a[k + 1 + k * a_dim1] != 0.)
					{
						k1 = k;
						k2 = k + 1;
						knext = k + 2;
					}
					else
					{
						k1 = k;
						k2 = k;
						knext = k + 1;
					}
				}

				if (l1 == l2 && k1 == k2)
				{
					i__3 = k1 - 1;
					suml = ddot_ (&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
					i__3 = l1 - 1;
					sumr = ddot_ (&i__3, &c__[k1 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
					vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
					scaloc = 1.;

					a11 = a[k1 + k1 * a_dim1] + sgn * b[l1 + l1 * b_dim1];
					da11 = abs (a11);
					if (da11 <= smin)
					{
						a11 = smin;
						da11 = smin;
						*info = 1;
					}
					db = abs (vec[0]);
					if (da11 < 1. && db > 1.)
					{
						if (db > bignum * da11)
						{
							scaloc = 1. / db;
						}
					}
					x[0] = vec[0] * scaloc / a11;

					if (scaloc != 1.)
					{
						i__3 = *n;
						for (j = 1; j <= i__3; ++j)
						{
							dscal_ (m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
							/* L70: */
						}
						*scale *= scaloc;
					}
					c__[k1 + l1 * c_dim1] = x[0];

				}
				else if (l1 == l2 && k1 != k2)
				{

					i__3 = k1 - 1;
					suml = ddot_ (&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
					i__3 = l1 - 1;
					sumr = ddot_ (&i__3, &c__[k1 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
					vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

					i__3 = k1 - 1;
					suml = ddot_ (&i__3, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
					i__3 = l1 - 1;
					sumr = ddot_ (&i__3, &c__[k2 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
					vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

					d__1 = -sgn * b[l1 + l1 * b_dim1];
					dlaln2_ (&c_true, &c__2, &c__1, &smin, &c_b348, &a[k1 + k1
																						* a_dim1], lda, &c_b348, &c_b348, vec, &c__2, &d__1, &c_b507, x, &c__2, &scaloc, &xnorm,
								&ierr);
					if (ierr != 0)
					{
						*info = 1;
					}

					if (scaloc != 1.)
					{
						i__3 = *n;
						for (j = 1; j <= i__3; ++j)
						{
							dscal_ (m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
							/* L80: */
						}
						*scale *= scaloc;
					}
					c__[k1 + l1 * c_dim1] = x[0];
					c__[k2 + l1 * c_dim1] = x[1];

				}
				else if (l1 != l2 && k1 == k2)
				{

					i__3 = k1 - 1;
					suml = ddot_ (&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
					i__3 = l1 - 1;
					sumr = ddot_ (&i__3, &c__[k1 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
					vec[0] = sgn * (c__[k1 + l1 * c_dim1] - (suml + sgn * sumr));

					i__3 = k1 - 1;
					suml = ddot_ (&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * c_dim1 + 1], &c__1);
					i__3 = l1 - 1;
					sumr = ddot_ (&i__3, &c__[k1 + c_dim1], ldc, &b[l2 * b_dim1 + 1], &c__1);
					vec[1] = sgn * (c__[k1 + l2 * c_dim1] - (suml + sgn * sumr));

					d__1 = -sgn * a[k1 + k1 * a_dim1];
					dlaln2_ (&c_true, &c__2, &c__1, &smin, &c_b348, &b[l1 + l1
																						* b_dim1], ldb, &c_b348, &c_b348, vec, &c__2, &d__1, &c_b507, x, &c__2, &scaloc, &xnorm,
								&ierr);
					if (ierr != 0)
					{
						*info = 1;
					}

					if (scaloc != 1.)
					{
						i__3 = *n;
						for (j = 1; j <= i__3; ++j)
						{
							dscal_ (m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
							/* L90: */
						}
						*scale *= scaloc;
					}
					c__[k1 + l1 * c_dim1] = x[0];
					c__[k1 + l2 * c_dim1] = x[1];

				}
				else if (l1 != l2 && k1 != k2)
				{

					i__3 = k1 - 1;
					suml = ddot_ (&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
					i__3 = l1 - 1;
					sumr = ddot_ (&i__3, &c__[k1 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
					vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

					i__3 = k1 - 1;
					suml = ddot_ (&i__3, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * c_dim1 + 1], &c__1);
					i__3 = l1 - 1;
					sumr = ddot_ (&i__3, &c__[k1 + c_dim1], ldc, &b[l2 * b_dim1 + 1], &c__1);
					vec[2] = c__[k1 + l2 * c_dim1] - (suml + sgn * sumr);

					i__3 = k1 - 1;
					suml = ddot_ (&i__3, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
					i__3 = l1 - 1;
					sumr = ddot_ (&i__3, &c__[k2 + c_dim1], ldc, &b[l1 * b_dim1 + 1], &c__1);
					vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

					i__3 = k1 - 1;
					suml = ddot_ (&i__3, &a[k2 * a_dim1 + 1], &c__1, &c__[l2 * c_dim1 + 1], &c__1);
					i__3 = l1 - 1;
					sumr = ddot_ (&i__3, &c__[k2 + c_dim1], ldc, &b[l2 * b_dim1 + 1], &c__1);
					vec[3] = c__[k2 + l2 * c_dim1] - (suml + sgn * sumr);

					dlasy2_ (&c_true, &c_false, isgn, &c__2, &c__2, &a[k1 + k1
																						* a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, vec, &c__2, &scaloc, x, &c__2, &xnorm, &ierr);
					if (ierr != 0)
					{
						*info = 1;
					}

					if (scaloc != 1.)
					{
						i__3 = *n;
						for (j = 1; j <= i__3; ++j)
						{
							dscal_ (m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
							/* L100: */
						}
						*scale *= scaloc;
					}
					c__[k1 + l1 * c_dim1] = x[0];
					c__[k1 + l2 * c_dim1] = x[2];
					c__[k2 + l1 * c_dim1] = x[1];
					c__[k2 + l2 * c_dim1] = x[3];
				}

			 L110:
				;
			}
		 L120:
			;
		}

	}
	else if (!notrna && !notrnb)
	{

		/*        Solve    A'*X + ISGN*X*B' = scale*C. */

		/*        The (K,L)th block of X is determined starting from */
		/*        top-right corner column by column by */

		/*           A(K,K)'*X(K,L) + ISGN*X(K,L)*B(L,L)' = C(K,L) - R(K,L) */

		/*        Where */
		/*                     K-1                          N */
		/*            R(K,L) = SUM [A(I,K)'*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)']. */
		/*                     I=1                        J=L+1 */

		/*        Start column loop (index = L) */
		/*        L1 (L2): column index of the first (last) row of X(K,L) */

		lnext = *n;
		for (l = *n; l >= 1; --l)
		{
			if (l > lnext)
			{
				goto L180;
			}
			if (l == 1)
			{
				l1 = l;
				l2 = l;
			}
			else
			{
				if (b[l + (l - 1) * b_dim1] != 0.)
				{
					l1 = l - 1;
					l2 = l;
					lnext = l - 2;
				}
				else
				{
					l1 = l;
					l2 = l;
					lnext = l - 1;
				}
			}

			/*           Start row loop (index = K) */
			/*           K1 (K2): row index of the first (last) row of X(K,L) */

			knext = 1;
			i__1 = *m;
			for (k = 1; k <= i__1; ++k)
			{
				if (k < knext)
				{
					goto L170;
				}
				if (k == *m)
				{
					k1 = k;
					k2 = k;
				}
				else
				{
					if (a[k + 1 + k * a_dim1] != 0.)
					{
						k1 = k;
						k2 = k + 1;
						knext = k + 2;
					}
					else
					{
						k1 = k;
						k2 = k;
						knext = k + 1;
					}
				}

				if (l1 == l2 && k1 == k2)
				{
					i__2 = k1 - 1;
					suml = ddot_ (&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
					i__2 = *n - l1;
					/* Computing MIN */
					i__3 = l1 + 1;
					/* Computing MIN */
					i__4 = l1 + 1;
					sumr = ddot_ (&i__2, &c__[k1 + min (i__3, *n) * c_dim1], ldc, &b[l1 + min (i__4, *n) * b_dim1], ldb);
					vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
					scaloc = 1.;

					a11 = a[k1 + k1 * a_dim1] + sgn * b[l1 + l1 * b_dim1];
					da11 = abs (a11);
					if (da11 <= smin)
					{
						a11 = smin;
						da11 = smin;
						*info = 1;
					}
					db = abs (vec[0]);
					if (da11 < 1. && db > 1.)
					{
						if (db > bignum * da11)
						{
							scaloc = 1. / db;
						}
					}
					x[0] = vec[0] * scaloc / a11;

					if (scaloc != 1.)
					{
						i__2 = *n;
						for (j = 1; j <= i__2; ++j)
						{
							dscal_ (m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
							/* L130: */
						}
						*scale *= scaloc;
					}
					c__[k1 + l1 * c_dim1] = x[0];

				}
				else if (l1 == l2 && k1 != k2)
				{

					i__2 = k1 - 1;
					suml = ddot_ (&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
					i__2 = *n - l2;
					/* Computing MIN */
					i__3 = l2 + 1;
					/* Computing MIN */
					i__4 = l2 + 1;
					sumr = ddot_ (&i__2, &c__[k1 + min (i__3, *n) * c_dim1], ldc, &b[l1 + min (i__4, *n) * b_dim1], ldb);
					vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

					i__2 = k1 - 1;
					suml = ddot_ (&i__2, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
					i__2 = *n - l2;
					/* Computing MIN */
					i__3 = l2 + 1;
					/* Computing MIN */
					i__4 = l2 + 1;
					sumr = ddot_ (&i__2, &c__[k2 + min (i__3, *n) * c_dim1], ldc, &b[l1 + min (i__4, *n) * b_dim1], ldb);
					vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

					d__1 = -sgn * b[l1 + l1 * b_dim1];
					dlaln2_ (&c_true, &c__2, &c__1, &smin, &c_b348, &a[k1 + k1
																						* a_dim1], lda, &c_b348, &c_b348, vec, &c__2, &d__1, &c_b507, x, &c__2, &scaloc, &xnorm,
								&ierr);
					if (ierr != 0)
					{
						*info = 1;
					}

					if (scaloc != 1.)
					{
						i__2 = *n;
						for (j = 1; j <= i__2; ++j)
						{
							dscal_ (m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
							/* L140: */
						}
						*scale *= scaloc;
					}
					c__[k1 + l1 * c_dim1] = x[0];
					c__[k2 + l1 * c_dim1] = x[1];

				}
				else if (l1 != l2 && k1 == k2)
				{

					i__2 = k1 - 1;
					suml = ddot_ (&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
					i__2 = *n - l2;
					/* Computing MIN */
					i__3 = l2 + 1;
					/* Computing MIN */
					i__4 = l2 + 1;
					sumr = ddot_ (&i__2, &c__[k1 + min (i__3, *n) * c_dim1], ldc, &b[l1 + min (i__4, *n) * b_dim1], ldb);
					vec[0] = sgn * (c__[k1 + l1 * c_dim1] - (suml + sgn * sumr));

					i__2 = k1 - 1;
					suml = ddot_ (&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * c_dim1 + 1], &c__1);
					i__2 = *n - l2;
					/* Computing MIN */
					i__3 = l2 + 1;
					/* Computing MIN */
					i__4 = l2 + 1;
					sumr = ddot_ (&i__2, &c__[k1 + min (i__3, *n) * c_dim1], ldc, &b[l2 + min (i__4, *n) * b_dim1], ldb);
					vec[1] = sgn * (c__[k1 + l2 * c_dim1] - (suml + sgn * sumr));

					d__1 = -sgn * a[k1 + k1 * a_dim1];
					dlaln2_ (&c_false, &c__2, &c__1, &smin, &c_b348, &b[l1 +
																						 l1 * b_dim1], ldb, &c_b348, &c_b348, vec, &c__2, &d__1, &c_b507, x, &c__2, &scaloc, &xnorm,
								&ierr);
					if (ierr != 0)
					{
						*info = 1;
					}

					if (scaloc != 1.)
					{
						i__2 = *n;
						for (j = 1; j <= i__2; ++j)
						{
							dscal_ (m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
							/* L150: */
						}
						*scale *= scaloc;
					}
					c__[k1 + l1 * c_dim1] = x[0];
					c__[k1 + l2 * c_dim1] = x[1];

				}
				else if (l1 != l2 && k1 != k2)
				{

					i__2 = k1 - 1;
					suml = ddot_ (&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
					i__2 = *n - l2;
					/* Computing MIN */
					i__3 = l2 + 1;
					/* Computing MIN */
					i__4 = l2 + 1;
					sumr = ddot_ (&i__2, &c__[k1 + min (i__3, *n) * c_dim1], ldc, &b[l1 + min (i__4, *n) * b_dim1], ldb);
					vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

					i__2 = k1 - 1;
					suml = ddot_ (&i__2, &a[k1 * a_dim1 + 1], &c__1, &c__[l2 * c_dim1 + 1], &c__1);
					i__2 = *n - l2;
					/* Computing MIN */
					i__3 = l2 + 1;
					/* Computing MIN */
					i__4 = l2 + 1;
					sumr = ddot_ (&i__2, &c__[k1 + min (i__3, *n) * c_dim1], ldc, &b[l2 + min (i__4, *n) * b_dim1], ldb);
					vec[2] = c__[k1 + l2 * c_dim1] - (suml + sgn * sumr);

					i__2 = k1 - 1;
					suml = ddot_ (&i__2, &a[k2 * a_dim1 + 1], &c__1, &c__[l1 * c_dim1 + 1], &c__1);
					i__2 = *n - l2;
					/* Computing MIN */
					i__3 = l2 + 1;
					/* Computing MIN */
					i__4 = l2 + 1;
					sumr = ddot_ (&i__2, &c__[k2 + min (i__3, *n) * c_dim1], ldc, &b[l1 + min (i__4, *n) * b_dim1], ldb);
					vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

					i__2 = k1 - 1;
					suml = ddot_ (&i__2, &a[k2 * a_dim1 + 1], &c__1, &c__[l2 * c_dim1 + 1], &c__1);
					i__2 = *n - l2;
					/* Computing MIN */
					i__3 = l2 + 1;
					/* Computing MIN */
					i__4 = l2 + 1;
					sumr = ddot_ (&i__2, &c__[k2 + min (i__3, *n) * c_dim1], ldc, &b[l2 + min (i__4, *n) * b_dim1], ldb);
					vec[3] = c__[k2 + l2 * c_dim1] - (suml + sgn * sumr);

					dlasy2_ (&c_true, &c_true, isgn, &c__2, &c__2, &a[k1 + k1 *
																					  a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, vec, &c__2, &scaloc, x, &c__2, &xnorm, &ierr);
					if (ierr != 0)
					{
						*info = 1;
					}

					if (scaloc != 1.)
					{
						i__2 = *n;
						for (j = 1; j <= i__2; ++j)
						{
							dscal_ (m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
							/* L160: */
						}
						*scale *= scaloc;
					}
					c__[k1 + l1 * c_dim1] = x[0];
					c__[k1 + l2 * c_dim1] = x[2];
					c__[k2 + l1 * c_dim1] = x[1];
					c__[k2 + l2 * c_dim1] = x[3];
				}

			 L170:
				;
			}
		 L180:
			;
		}

	}
	else if (notrna && !notrnb)
	{

		/*        Solve    A*X + ISGN*X*B' = scale*C. */

		/*        The (K,L)th block of X is determined starting from */
		/*        bottom-right corner column by column by */

		/*            A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)' = C(K,L) - R(K,L) */

		/*        Where */
		/*                      M                          N */
		/*            R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)']. */
		/*                    I=K+1                      J=L+1 */

		/*        Start column loop (index = L) */
		/*        L1 (L2): column index of the first (last) row of X(K,L) */

		lnext = *n;
		for (l = *n; l >= 1; --l)
		{
			if (l > lnext)
			{
				goto L240;
			}
			if (l == 1)
			{
				l1 = l;
				l2 = l;
			}
			else
			{
				if (b[l + (l - 1) * b_dim1] != 0.)
				{
					l1 = l - 1;
					l2 = l;
					lnext = l - 2;
				}
				else
				{
					l1 = l;
					l2 = l;
					lnext = l - 1;
				}
			}

			/*           Start row loop (index = K) */
			/*           K1 (K2): row index of the first (last) row of X(K,L) */

			knext = *m;
			for (k = *m; k >= 1; --k)
			{
				if (k > knext)
				{
					goto L230;
				}
				if (k == 1)
				{
					k1 = k;
					k2 = k;
				}
				else
				{
					if (a[k + (k - 1) * a_dim1] != 0.)
					{
						k1 = k - 1;
						k2 = k;
						knext = k - 2;
					}
					else
					{
						k1 = k;
						k2 = k;
						knext = k - 1;
					}
				}

				if (l1 == l2 && k1 == k2)
				{
					i__1 = *m - k1;
					/* Computing MIN */
					i__2 = k1 + 1;
					/* Computing MIN */
					i__3 = k1 + 1;
					suml = ddot_ (&i__1, &a[k1 + min (i__2, *m) * a_dim1], lda, &c__[min (i__3, *m) + l1 * c_dim1], &c__1);
					i__1 = *n - l1;
					/* Computing MIN */
					i__2 = l1 + 1;
					/* Computing MIN */
					i__3 = l1 + 1;
					sumr = ddot_ (&i__1, &c__[k1 + min (i__2, *n) * c_dim1], ldc, &b[l1 + min (i__3, *n) * b_dim1], ldb);
					vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);
					scaloc = 1.;

					a11 = a[k1 + k1 * a_dim1] + sgn * b[l1 + l1 * b_dim1];
					da11 = abs (a11);
					if (da11 <= smin)
					{
						a11 = smin;
						da11 = smin;
						*info = 1;
					}
					db = abs (vec[0]);
					if (da11 < 1. && db > 1.)
					{
						if (db > bignum * da11)
						{
							scaloc = 1. / db;
						}
					}
					x[0] = vec[0] * scaloc / a11;

					if (scaloc != 1.)
					{
						i__1 = *n;
						for (j = 1; j <= i__1; ++j)
						{
							dscal_ (m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
							/* L190: */
						}
						*scale *= scaloc;
					}
					c__[k1 + l1 * c_dim1] = x[0];

				}
				else if (l1 == l2 && k1 != k2)
				{

					i__1 = *m - k2;
					/* Computing MIN */
					i__2 = k2 + 1;
					/* Computing MIN */
					i__3 = k2 + 1;
					suml = ddot_ (&i__1, &a[k1 + min (i__2, *m) * a_dim1], lda, &c__[min (i__3, *m) + l1 * c_dim1], &c__1);
					i__1 = *n - l2;
					/* Computing MIN */
					i__2 = l2 + 1;
					/* Computing MIN */
					i__3 = l2 + 1;
					sumr = ddot_ (&i__1, &c__[k1 + min (i__2, *n) * c_dim1], ldc, &b[l1 + min (i__3, *n) * b_dim1], ldb);
					vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

					i__1 = *m - k2;
					/* Computing MIN */
					i__2 = k2 + 1;
					/* Computing MIN */
					i__3 = k2 + 1;
					suml = ddot_ (&i__1, &a[k2 + min (i__2, *m) * a_dim1], lda, &c__[min (i__3, *m) + l1 * c_dim1], &c__1);
					i__1 = *n - l2;
					/* Computing MIN */
					i__2 = l2 + 1;
					/* Computing MIN */
					i__3 = l2 + 1;
					sumr = ddot_ (&i__1, &c__[k2 + min (i__2, *n) * c_dim1], ldc, &b[l1 + min (i__3, *n) * b_dim1], ldb);
					vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

					d__1 = -sgn * b[l1 + l1 * b_dim1];
					dlaln2_ (&c_false, &c__2, &c__1, &smin, &c_b348, &a[k1 +
																						 k1 * a_dim1], lda, &c_b348, &c_b348, vec, &c__2, &d__1, &c_b507, x, &c__2, &scaloc, &xnorm,
								&ierr);
					if (ierr != 0)
					{
						*info = 1;
					}

					if (scaloc != 1.)
					{
						i__1 = *n;
						for (j = 1; j <= i__1; ++j)
						{
							dscal_ (m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
							/* L200: */
						}
						*scale *= scaloc;
					}
					c__[k1 + l1 * c_dim1] = x[0];
					c__[k2 + l1 * c_dim1] = x[1];

				}
				else if (l1 != l2 && k1 == k2)
				{

					i__1 = *m - k1;
					/* Computing MIN */
					i__2 = k1 + 1;
					/* Computing MIN */
					i__3 = k1 + 1;
					suml = ddot_ (&i__1, &a[k1 + min (i__2, *m) * a_dim1], lda, &c__[min (i__3, *m) + l1 * c_dim1], &c__1);
					i__1 = *n - l2;
					/* Computing MIN */
					i__2 = l2 + 1;
					/* Computing MIN */
					i__3 = l2 + 1;
					sumr = ddot_ (&i__1, &c__[k1 + min (i__2, *n) * c_dim1], ldc, &b[l1 + min (i__3, *n) * b_dim1], ldb);
					vec[0] = sgn * (c__[k1 + l1 * c_dim1] - (suml + sgn * sumr));

					i__1 = *m - k1;
					/* Computing MIN */
					i__2 = k1 + 1;
					/* Computing MIN */
					i__3 = k1 + 1;
					suml = ddot_ (&i__1, &a[k1 + min (i__2, *m) * a_dim1], lda, &c__[min (i__3, *m) + l2 * c_dim1], &c__1);
					i__1 = *n - l2;
					/* Computing MIN */
					i__2 = l2 + 1;
					/* Computing MIN */
					i__3 = l2 + 1;
					sumr = ddot_ (&i__1, &c__[k1 + min (i__2, *n) * c_dim1], ldc, &b[l2 + min (i__3, *n) * b_dim1], ldb);
					vec[1] = sgn * (c__[k1 + l2 * c_dim1] - (suml + sgn * sumr));

					d__1 = -sgn * a[k1 + k1 * a_dim1];
					dlaln2_ (&c_false, &c__2, &c__1, &smin, &c_b348, &b[l1 +
																						 l1 * b_dim1], ldb, &c_b348, &c_b348, vec, &c__2, &d__1, &c_b507, x, &c__2, &scaloc, &xnorm,
								&ierr);
					if (ierr != 0)
					{
						*info = 1;
					}

					if (scaloc != 1.)
					{
						i__1 = *n;
						for (j = 1; j <= i__1; ++j)
						{
							dscal_ (m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
							/* L210: */
						}
						*scale *= scaloc;
					}
					c__[k1 + l1 * c_dim1] = x[0];
					c__[k1 + l2 * c_dim1] = x[1];

				}
				else if (l1 != l2 && k1 != k2)
				{

					i__1 = *m - k2;
					/* Computing MIN */
					i__2 = k2 + 1;
					/* Computing MIN */
					i__3 = k2 + 1;
					suml = ddot_ (&i__1, &a[k1 + min (i__2, *m) * a_dim1], lda, &c__[min (i__3, *m) + l1 * c_dim1], &c__1);
					i__1 = *n - l2;
					/* Computing MIN */
					i__2 = l2 + 1;
					/* Computing MIN */
					i__3 = l2 + 1;
					sumr = ddot_ (&i__1, &c__[k1 + min (i__2, *n) * c_dim1], ldc, &b[l1 + min (i__3, *n) * b_dim1], ldb);
					vec[0] = c__[k1 + l1 * c_dim1] - (suml + sgn * sumr);

					i__1 = *m - k2;
					/* Computing MIN */
					i__2 = k2 + 1;
					/* Computing MIN */
					i__3 = k2 + 1;
					suml = ddot_ (&i__1, &a[k1 + min (i__2, *m) * a_dim1], lda, &c__[min (i__3, *m) + l2 * c_dim1], &c__1);
					i__1 = *n - l2;
					/* Computing MIN */
					i__2 = l2 + 1;
					/* Computing MIN */
					i__3 = l2 + 1;
					sumr = ddot_ (&i__1, &c__[k1 + min (i__2, *n) * c_dim1], ldc, &b[l2 + min (i__3, *n) * b_dim1], ldb);
					vec[2] = c__[k1 + l2 * c_dim1] - (suml + sgn * sumr);

					i__1 = *m - k2;
					/* Computing MIN */
					i__2 = k2 + 1;
					/* Computing MIN */
					i__3 = k2 + 1;
					suml = ddot_ (&i__1, &a[k2 + min (i__2, *m) * a_dim1], lda, &c__[min (i__3, *m) + l1 * c_dim1], &c__1);
					i__1 = *n - l2;
					/* Computing MIN */
					i__2 = l2 + 1;
					/* Computing MIN */
					i__3 = l2 + 1;
					sumr = ddot_ (&i__1, &c__[k2 + min (i__2, *n) * c_dim1], ldc, &b[l1 + min (i__3, *n) * b_dim1], ldb);
					vec[1] = c__[k2 + l1 * c_dim1] - (suml + sgn * sumr);

					i__1 = *m - k2;
					/* Computing MIN */
					i__2 = k2 + 1;
					/* Computing MIN */
					i__3 = k2 + 1;
					suml = ddot_ (&i__1, &a[k2 + min (i__2, *m) * a_dim1], lda, &c__[min (i__3, *m) + l2 * c_dim1], &c__1);
					i__1 = *n - l2;
					/* Computing MIN */
					i__2 = l2 + 1;
					/* Computing MIN */
					i__3 = l2 + 1;
					sumr = ddot_ (&i__1, &c__[k2 + min (i__2, *n) * c_dim1], ldc, &b[l2 + min (i__3, *n) * b_dim1], ldb);
					vec[3] = c__[k2 + l2 * c_dim1] - (suml + sgn * sumr);

					dlasy2_ (&c_false, &c_true, isgn, &c__2, &c__2, &a[k1 + k1
																						* a_dim1], lda, &b[l1 + l1 * b_dim1], ldb, vec, &c__2, &scaloc, x, &c__2, &xnorm, &ierr);
					if (ierr != 0)
					{
						*info = 1;
					}

					if (scaloc != 1.)
					{
						i__1 = *n;
						for (j = 1; j <= i__1; ++j)
						{
							dscal_ (m, &scaloc, &c__[j * c_dim1 + 1], &c__1);
							/* L220: */
						}
						*scale *= scaloc;
					}
					c__[k1 + l1 * c_dim1] = x[0];
					c__[k1 + l2 * c_dim1] = x[2];
					c__[k2 + l1 * c_dim1] = x[1];
					c__[k2 + l2 * c_dim1] = x[3];
				}

			 L230:
				;
			}
		 L240:
			;
		}

	}

	return 0;

	/*     End of DTRSYL */

}	/* dtrsyl_ */

static integer
ieeeck_ (integer * ispec, real * zero, real * one)
{
	/* System generated locals */
	integer ret_val;

	/* Local variables */
	static real nan1, nan2, nan3, nan4, nan5, nan6, neginf, posinf, negzro, newzro;


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     June 30, 1998 */

	/*     .. Scalar Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  IEEECK is called from the ILAENV to verify that Infinity and */
	/*  possibly NaN arithmetic is safe (i.e. will not trap). */

	/*  Arguments */
	/*  ========= */

	/*  ISPEC   (input) INTEGER */
	/*          Specifies whether to test just for inifinity arithmetic */
	/*          or whether to test for infinity and NaN arithmetic. */
	/*          = 0: Verify infinity arithmetic only. */
	/*          = 1: Verify infinity and NaN arithmetic. */

	/*  ZERO    (input) REAL */
	/*          Must contain the value 0.0 */
	/*          This is passed to prevent the compiler from optimizing */
	/*          away this code. */

	/*  ONE     (input) REAL */
	/*          Must contain the value 1.0 */
	/*          This is passed to prevent the compiler from optimizing */
	/*          away this code. */

	/*  RETURN VALUE:  INTEGER */
	/*          = 0:  Arithmetic failed to produce the correct answers */
	/*          = 1:  Arithmetic produced the correct answers */

	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Executable Statements .. */
	ret_val = 1;

	posinf = *one / *zero;
	if (posinf <= *one)
	{
		ret_val = 0;
		return ret_val;
	}

	neginf = -(*one) / *zero;
	if (neginf >= *zero)
	{
		ret_val = 0;
		return ret_val;
	}

	negzro = *one / (neginf + *one);
	if (negzro != *zero)
	{
		ret_val = 0;
		return ret_val;
	}

	neginf = *one / negzro;
	if (neginf >= *zero)
	{
		ret_val = 0;
		return ret_val;
	}

	newzro = negzro + *zero;
	if (newzro != *zero)
	{
		ret_val = 0;
		return ret_val;
	}

	posinf = *one / newzro;
	if (posinf <= *one)
	{
		ret_val = 0;
		return ret_val;
	}

	neginf *= posinf;
	if (neginf >= *zero)
	{
		ret_val = 0;
		return ret_val;
	}

	posinf *= posinf;
	if (posinf <= *one)
	{
		ret_val = 0;
		return ret_val;
	}




	/*     Return if we were only asked to check infinity arithmetic */

	if (*ispec == 0)
	{
		return ret_val;
	}

	nan1 = posinf + neginf;

	nan2 = posinf / neginf;

	nan3 = posinf / posinf;

	nan4 = posinf * *zero;

	nan5 = neginf * negzro;

	nan6 = nan5 * 0.f;

	if (nan1 == nan1)
	{
		ret_val = 0;
		return ret_val;
	}

	if (nan2 == nan2)
	{
		ret_val = 0;
		return ret_val;
	}

	if (nan3 == nan3)
	{
		ret_val = 0;
		return ret_val;
	}

	if (nan4 == nan4)
	{
		ret_val = 0;
		return ret_val;
	}

	if (nan5 == nan5)
	{
		ret_val = 0;
		return ret_val;
	}

	if (nan6 == nan6)
	{
		ret_val = 0;
		return ret_val;
	}

	return ret_val;
}	/* ieeeck_ */

static integer
ilaenv_ (integer * ispec, char *name__, char *opts, integer * n1, integer * n2, integer * n3, integer * n4, ftnlen name_len, ftnlen opts_len)
{
	/* System generated locals */
	integer ret_val;

	/* Builtin functions */
	/* Subroutine */ int s_copy (char *, char *, ftnlen, ftnlen);
	integer s_cmp (char *, char *, ftnlen, ftnlen);

	/* Local variables */
	static integer i__;
	static char c1[1], c2[2], c3[3], c4[2];
	static integer ic, nb, iz, nx;
	static logical cname, sname;
	static integer nbmin;
	static char subnam[6];


	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     June 30, 1999 */

	/*     .. Scalar Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  ILAENV is called from the LAPACK routines to choose problem-dependent */
	/*  parameters for the local environment.  See ISPEC for a description of */
	/*  the parameters. */

	/*  This version provides a set of parameters which should give good, */
	/*  but not optimal, performance on many of the currently available */
	/*  computers.  Users are encouraged to modify this subroutine to set */
	/*  the tuning parameters for their particular machine using the option */
	/*  and problem size information in the arguments. */

	/*  This routine will not function correctly if it is converted to all */
	/*  lower case.  Converting it to all upper case is allowed. */

	/*  Arguments */
	/*  ========= */

	/*  ISPEC   (input) INTEGER */
	/*          Specifies the parameter to be returned as the value of */
	/*          ILAENV. */
	/*          = 1: the optimal blocksize; if this value is 1, an unblocked */
	/*               algorithm will give the best performance. */
	/*          = 2: the minimum block size for which the block routine */
	/*               should be used; if the usable block size is less than */
	/*               this value, an unblocked routine should be used. */
	/*          = 3: the crossover point (in a block routine, for N less */
	/*               than this value, an unblocked routine should be used) */
	/*          = 4: the number of shifts, used in the nonsymmetric */
	/*               eigenvalue routines */
	/*          = 5: the minimum column dimension for blocking to be used; */
	/*               rectangular blocks must have dimension at least k by m, */
	/*               where k is given by ILAENV(2,...) and m by ILAENV(5,...) */
	/*          = 6: the crossover point for the SVD (when reducing an m by n */
	/*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds */
	/*               this value, a QR factorization is used first to reduce */
	/*               the matrix to a triangular form.) */
	/*          = 7: the number of processors */
	/*          = 8: the crossover point for the multishift QR and QZ methods */
	/*               for nonsymmetric eigenvalue problems. */
	/*          = 9: maximum size of the subproblems at the bottom of the */
	/*               computation tree in the divide-and-conquer algorithm */
	/*               (used by xGELSD and xGESDD) */
	/*          =10: ieee NaN arithmetic can be trusted not to trap */
	/*          =11: infinity arithmetic can be trusted not to trap */

	/*  NAME    (input) CHARACTER*(*) */
	/*          The name of the calling subroutine, in either upper case or */
	/*          lower case. */

	/*  OPTS    (input) CHARACTER*(*) */
	/*          The character options to the subroutine NAME, concatenated */
	/*          into a single character string.  For example, UPLO = 'U', */
	/*          TRANS = 'T', and DIAG = 'N' for a triangular routine would */
	/*          be specified as OPTS = 'UTN'. */

	/*  N1      (input) INTEGER */
	/*  N2      (input) INTEGER */
	/*  N3      (input) INTEGER */
	/*  N4      (input) INTEGER */
	/*          Problem dimensions for the subroutine NAME; these may not all */
	/*          be required. */

	/* (ILAENV) (output) INTEGER */
	/*          >= 0: the value of the parameter specified by ISPEC */
	/*          < 0:  if ILAENV = -k, the k-th argument had an illegal value. */

	/*  Further Details */
	/*  =============== */

	/*  The following conventions have been used when calling ILAENV from the */
	/*  LAPACK routines: */
	/*  1)  OPTS is a concatenation of all of the character options to */
	/*      subroutine NAME, in the same order that they appear in the */
	/*      argument list for NAME, even if they are not used in determining */
	/*      the value of the parameter specified by ISPEC. */
	/*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order */
	/*      that they appear in the argument list for NAME.  N1 is used */
	/*      first, N2 second, and so on, and unused problem dimensions are */
	/*      passed a value of -1. */
	/*  3)  The parameter value returned by ILAENV is checked for validity in */
	/*      the calling subroutine.  For example, ILAENV is used to retrieve */
	/*      the optimal blocksize for STRTRI as follows: */

	/*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 ) */
	/*      IF( NB.LE.1 ) NB = MAX( 1, N ) */

	/*  ===================================================================== */

	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	switch (*ispec)
	{
	case 1:
		goto L100;
	case 2:
		goto L100;
	case 3:
		goto L100;
	case 4:
		goto L400;
	case 5:
		goto L500;
	case 6:
		goto L600;
	case 7:
		goto L700;
	case 8:
		goto L800;
	case 9:
		goto L900;
	case 10:
		goto L1000;
	case 11:
		goto L1100;
	}

	/*     Invalid value for ISPEC */

	ret_val = -1;
	return ret_val;

 L100:

	/*     Convert NAME to upper case if the first character is lower case. */

	ret_val = 1;
	s_copy (subnam, name__, (ftnlen) 6, name_len);
	ic = *(unsigned char *) subnam;
	iz = 'Z';
	if (iz == 90 || iz == 122)
	{

		/*        ASCII character set */

		if (ic >= 97 && ic <= 122)
		{
			*(unsigned char *) subnam = (char) (ic - 32);
			for (i__ = 2; i__ <= 6; ++i__)
			{
				ic = *(unsigned char *) &subnam[i__ - 1];
				if (ic >= 97 && ic <= 122)
				{
					*(unsigned char *) &subnam[i__ - 1] = (char) (ic - 32);
				}
				/* L10: */
			}
		}

	}
	else if (iz == 233 || iz == 169)
	{

		/*        EBCDIC character set */

		if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 && ic <= 169)
		{
			*(unsigned char *) subnam = (char) (ic + 64);
			for (i__ = 2; i__ <= 6; ++i__)
			{
				ic = *(unsigned char *) &subnam[i__ - 1];
				if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 && ic <= 169)
				{
					*(unsigned char *) &subnam[i__ - 1] = (char) (ic + 64);
				}
				/* L20: */
			}
		}

	}
	else if (iz == 218 || iz == 250)
	{

		/*        Prime machines:  ASCII+128 */

		if (ic >= 225 && ic <= 250)
		{
			*(unsigned char *) subnam = (char) (ic - 32);
			for (i__ = 2; i__ <= 6; ++i__)
			{
				ic = *(unsigned char *) &subnam[i__ - 1];
				if (ic >= 225 && ic <= 250)
				{
					*(unsigned char *) &subnam[i__ - 1] = (char) (ic - 32);
				}
				/* L30: */
			}
		}
	}

	*(unsigned char *) c1 = *(unsigned char *) subnam;
	sname = *(unsigned char *) c1 == 'S' || *(unsigned char *) c1 == 'D';
	cname = *(unsigned char *) c1 == 'C' || *(unsigned char *) c1 == 'Z';
	if (!(cname || sname))
	{
		return ret_val;
	}
	s_copy (c2, subnam + 1, (ftnlen) 2, (ftnlen) 2);
	s_copy (c3, subnam + 3, (ftnlen) 3, (ftnlen) 3);
	s_copy (c4, c3 + 1, (ftnlen) 2, (ftnlen) 2);

	switch (*ispec)
	{
	case 1:
		goto L110;
	case 2:
		goto L200;
	case 3:
		goto L300;
	}

 L110:

	/*     ISPEC = 1:  block size */

	/*     In these examples, separate code is provided for setting NB for */
	/*     real and complex.  We assume that NB will take the same value in */
	/*     single or double precision. */

	nb = 1;

	if (s_cmp (c2, "GE", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (s_cmp (c3, "TRF", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			if (sname)
			{
				nb = 64;
			}
			else
			{
				nb = 64;
			}
		}
		else if (s_cmp (c3, "QRF", (ftnlen) 3, (ftnlen) 3) == 0 || s_cmp (c3,
																								"RQF", (ftnlen) 3, (ftnlen) 3) == 0 || s_cmp (c3, "LQF", (ftnlen)
																																							 3, (ftnlen) 3) == 0
					|| s_cmp (c3, "QLF", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			if (sname)
			{
				nb = 32;
			}
			else
			{
				nb = 32;
			}
		}
		else if (s_cmp (c3, "HRD", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			if (sname)
			{
				nb = 32;
			}
			else
			{
				nb = 32;
			}
		}
		else if (s_cmp (c3, "BRD", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			if (sname)
			{
				nb = 32;
			}
			else
			{
				nb = 32;
			}
		}
		else if (s_cmp (c3, "TRI", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			if (sname)
			{
				nb = 64;
			}
			else
			{
				nb = 64;
			}
		}
	}
	else if (s_cmp (c2, "PO", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (s_cmp (c3, "TRF", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			if (sname)
			{
				nb = 64;
			}
			else
			{
				nb = 64;
			}
		}
	}
	else if (s_cmp (c2, "SY", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (s_cmp (c3, "TRF", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			if (sname)
			{
				nb = 64;
			}
			else
			{
				nb = 64;
			}
		}
		else if (sname && s_cmp (c3, "TRD", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			nb = 32;
		}
		else if (sname && s_cmp (c3, "GST", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			nb = 64;
		}
	}
	else if (cname && s_cmp (c2, "HE", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (s_cmp (c3, "TRF", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			nb = 64;
		}
		else if (s_cmp (c3, "TRD", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			nb = 32;
		}
		else if (s_cmp (c3, "GST", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			nb = 64;
		}
	}
	else if (sname && s_cmp (c2, "OR", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (*(unsigned char *) c3 == 'G')
		{
			if (s_cmp (c4, "QR", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "RQ",
																							(ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "LQ", (ftnlen) 2, (ftnlen) 2) == 0
				 || s_cmp (c4, "QL", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "HR", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "TR", (ftnlen) 2, (ftnlen) 2) == 0
				 || s_cmp (c4, "BR", (ftnlen) 2, (ftnlen) 2) == 0)
			{
				nb = 32;
			}
		}
		else if (*(unsigned char *) c3 == 'M')
		{
			if (s_cmp (c4, "QR", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "RQ",
																							(ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "LQ", (ftnlen) 2, (ftnlen) 2) == 0
				 || s_cmp (c4, "QL", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "HR", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "TR", (ftnlen) 2, (ftnlen) 2) == 0
				 || s_cmp (c4, "BR", (ftnlen) 2, (ftnlen) 2) == 0)
			{
				nb = 32;
			}
		}
	}
	else if (cname && s_cmp (c2, "UN", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (*(unsigned char *) c3 == 'G')
		{
			if (s_cmp (c4, "QR", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "RQ",
																							(ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "LQ", (ftnlen) 2, (ftnlen) 2) == 0
				 || s_cmp (c4, "QL", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "HR", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "TR", (ftnlen) 2, (ftnlen) 2) == 0
				 || s_cmp (c4, "BR", (ftnlen) 2, (ftnlen) 2) == 0)
			{
				nb = 32;
			}
		}
		else if (*(unsigned char *) c3 == 'M')
		{
			if (s_cmp (c4, "QR", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "RQ",
																							(ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "LQ", (ftnlen) 2, (ftnlen) 2) == 0
				 || s_cmp (c4, "QL", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "HR", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "TR", (ftnlen) 2, (ftnlen) 2) == 0
				 || s_cmp (c4, "BR", (ftnlen) 2, (ftnlen) 2) == 0)
			{
				nb = 32;
			}
		}
	}
	else if (s_cmp (c2, "GB", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (s_cmp (c3, "TRF", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			if (sname)
			{
				if (*n4 <= 64)
				{
					nb = 1;
				}
				else
				{
					nb = 32;
				}
			}
			else
			{
				if (*n4 <= 64)
				{
					nb = 1;
				}
				else
				{
					nb = 32;
				}
			}
		}
	}
	else if (s_cmp (c2, "PB", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (s_cmp (c3, "TRF", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			if (sname)
			{
				if (*n2 <= 64)
				{
					nb = 1;
				}
				else
				{
					nb = 32;
				}
			}
			else
			{
				if (*n2 <= 64)
				{
					nb = 1;
				}
				else
				{
					nb = 32;
				}
			}
		}
	}
	else if (s_cmp (c2, "TR", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (s_cmp (c3, "TRI", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			if (sname)
			{
				nb = 64;
			}
			else
			{
				nb = 64;
			}
		}
	}
	else if (s_cmp (c2, "LA", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (s_cmp (c3, "UUM", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			if (sname)
			{
				nb = 64;
			}
			else
			{
				nb = 64;
			}
		}
	}
	else if (sname && s_cmp (c2, "ST", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (s_cmp (c3, "EBZ", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			nb = 1;
		}
	}
	ret_val = nb;
	return ret_val;

 L200:

	/*     ISPEC = 2:  minimum block size */

	nbmin = 2;
	if (s_cmp (c2, "GE", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (s_cmp (c3, "QRF", (ftnlen) 3, (ftnlen) 3) == 0 || s_cmp (c3, "RQF", (ftnlen) 3, (ftnlen) 3) == 0 || s_cmp (c3, "LQF", (ftnlen) 3, (ftnlen) 3) == 0
			 || s_cmp (c3, "QLF", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			if (sname)
			{
				nbmin = 2;
			}
			else
			{
				nbmin = 2;
			}
		}
		else if (s_cmp (c3, "HRD", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			if (sname)
			{
				nbmin = 2;
			}
			else
			{
				nbmin = 2;
			}
		}
		else if (s_cmp (c3, "BRD", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			if (sname)
			{
				nbmin = 2;
			}
			else
			{
				nbmin = 2;
			}
		}
		else if (s_cmp (c3, "TRI", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			if (sname)
			{
				nbmin = 2;
			}
			else
			{
				nbmin = 2;
			}
		}
	}
	else if (s_cmp (c2, "SY", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (s_cmp (c3, "TRF", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			if (sname)
			{
				nbmin = 8;
			}
			else
			{
				nbmin = 8;
			}
		}
		else if (sname && s_cmp (c3, "TRD", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			nbmin = 2;
		}
	}
	else if (cname && s_cmp (c2, "HE", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (s_cmp (c3, "TRD", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			nbmin = 2;
		}
	}
	else if (sname && s_cmp (c2, "OR", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (*(unsigned char *) c3 == 'G')
		{
			if (s_cmp (c4, "QR", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "RQ",
																							(ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "LQ", (ftnlen) 2, (ftnlen) 2) == 0
				 || s_cmp (c4, "QL", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "HR", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "TR", (ftnlen) 2, (ftnlen) 2) == 0
				 || s_cmp (c4, "BR", (ftnlen) 2, (ftnlen) 2) == 0)
			{
				nbmin = 2;
			}
		}
		else if (*(unsigned char *) c3 == 'M')
		{
			if (s_cmp (c4, "QR", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "RQ",
																							(ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "LQ", (ftnlen) 2, (ftnlen) 2) == 0
				 || s_cmp (c4, "QL", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "HR", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "TR", (ftnlen) 2, (ftnlen) 2) == 0
				 || s_cmp (c4, "BR", (ftnlen) 2, (ftnlen) 2) == 0)
			{
				nbmin = 2;
			}
		}
	}
	else if (cname && s_cmp (c2, "UN", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (*(unsigned char *) c3 == 'G')
		{
			if (s_cmp (c4, "QR", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "RQ",
																							(ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "LQ", (ftnlen) 2, (ftnlen) 2) == 0
				 || s_cmp (c4, "QL", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "HR", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "TR", (ftnlen) 2, (ftnlen) 2) == 0
				 || s_cmp (c4, "BR", (ftnlen) 2, (ftnlen) 2) == 0)
			{
				nbmin = 2;
			}
		}
		else if (*(unsigned char *) c3 == 'M')
		{
			if (s_cmp (c4, "QR", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "RQ",
																							(ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "LQ", (ftnlen) 2, (ftnlen) 2) == 0
				 || s_cmp (c4, "QL", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "HR", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "TR", (ftnlen) 2, (ftnlen) 2) == 0
				 || s_cmp (c4, "BR", (ftnlen) 2, (ftnlen) 2) == 0)
			{
				nbmin = 2;
			}
		}
	}
	ret_val = nbmin;
	return ret_val;

 L300:

	/*     ISPEC = 3:  crossover point */

	nx = 0;
	if (s_cmp (c2, "GE", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (s_cmp (c3, "QRF", (ftnlen) 3, (ftnlen) 3) == 0 || s_cmp (c3, "RQF", (ftnlen) 3, (ftnlen) 3) == 0 || s_cmp (c3, "LQF", (ftnlen) 3, (ftnlen) 3) == 0
			 || s_cmp (c3, "QLF", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			if (sname)
			{
				nx = 128;
			}
			else
			{
				nx = 128;
			}
		}
		else if (s_cmp (c3, "HRD", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			if (sname)
			{
				nx = 128;
			}
			else
			{
				nx = 128;
			}
		}
		else if (s_cmp (c3, "BRD", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			if (sname)
			{
				nx = 128;
			}
			else
			{
				nx = 128;
			}
		}
	}
	else if (s_cmp (c2, "SY", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (sname && s_cmp (c3, "TRD", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			nx = 32;
		}
	}
	else if (cname && s_cmp (c2, "HE", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (s_cmp (c3, "TRD", (ftnlen) 3, (ftnlen) 3) == 0)
		{
			nx = 32;
		}
	}
	else if (sname && s_cmp (c2, "OR", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (*(unsigned char *) c3 == 'G')
		{
			if (s_cmp (c4, "QR", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "RQ",
																							(ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "LQ", (ftnlen) 2, (ftnlen) 2) == 0
				 || s_cmp (c4, "QL", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "HR", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "TR", (ftnlen) 2, (ftnlen) 2) == 0
				 || s_cmp (c4, "BR", (ftnlen) 2, (ftnlen) 2) == 0)
			{
				nx = 128;
			}
		}
	}
	else if (cname && s_cmp (c2, "UN", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		if (*(unsigned char *) c3 == 'G')
		{
			if (s_cmp (c4, "QR", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "RQ",
																							(ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "LQ", (ftnlen) 2, (ftnlen) 2) == 0
				 || s_cmp (c4, "QL", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "HR", (ftnlen) 2, (ftnlen) 2) == 0 || s_cmp (c4, "TR", (ftnlen) 2, (ftnlen) 2) == 0
				 || s_cmp (c4, "BR", (ftnlen) 2, (ftnlen) 2) == 0)
			{
				nx = 128;
			}
		}
	}
	ret_val = nx;
	return ret_val;

 L400:

	/*     ISPEC = 4:  number of shifts (used by xHSEQR) */

	ret_val = 6;
	return ret_val;

 L500:

	/*     ISPEC = 5:  minimum column dimension (not used) */

	ret_val = 2;
	return ret_val;

 L600:

	/*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD) */

	ret_val = (integer) ((real) min (*n1, *n2) * 1.6f);
	return ret_val;

 L700:

	/*     ISPEC = 7:  number of processors (not used) */

	ret_val = 1;
	return ret_val;

 L800:

	/*     ISPEC = 8:  crossover point for multishift (used by xHSEQR) */

	ret_val = 50;
	return ret_val;

 L900:

	/*     ISPEC = 9:  maximum size of the subproblems at the bottom of the */
	/*                 computation tree in the divide-and-conquer algorithm */
	/*                 (used by xGELSD and xGESDD) */

	ret_val = 25;
	return ret_val;

 L1000:

	/*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap */

	/*     ILAENV = 0 */
	ret_val = 1;
	if (ret_val == 1)
	{
		ret_val = ieeeck_ (&c__0, &c_b2255, &c_b2256);
	}
	return ret_val;

 L1100:

	/*     ISPEC = 11: infinity arithmetic can be trusted not to trap */

	/*     ILAENV = 0 */
	ret_val = 1;
	if (ret_val == 1)
	{
		ret_val = ieeeck_ (&c__1, &c_b2255, &c_b2256);
	}
	return ret_val;

	/*     End of ILAENV */

}	/* ilaenv_ */

// logical lsame_(char *ca, char *cb, ftnlen ca_len, ftnlen cb_len)
// {
//    logical ret_val;
// 
//    static integer inta, intb, zcode;
// 
// 
//    ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
//    if (ret_val)
//    {
//       return ret_val;
//    }
// 
//    /*     Now test for equivalence if both characters are alphabetic. */
// 
//    zcode = 'Z';
// 
//    /*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime */
//    /*     machines, on which ICHAR returns a value with bit 8 set. */
//    /*     ICHAR('A') on Prime machines returns 193 which is the same as */
//    /*     ICHAR('A') on an EBCDIC machine. */
// 
//    inta = *(unsigned char *)ca;
//    intb = *(unsigned char *)cb;
// 
//    if (zcode == 90 || zcode == 122)
//    {
//       /*        ASCII is assumed - ZCODE is the ASCII code of either lower or */
//       /*        upper case 'Z'. */
//       if (inta >= 97 && inta <= 122) 
//       {
//          inta += -32;
//       }
//       if (intb >= 97 && intb <= 122) 
//       {
//          intb += -32;
//       }
//    } else if (zcode == 233 || zcode == 169) 
//    {
//       /*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or */
//       /*        upper case 'Z'. */
//       if (inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta >= 162 && inta <= 169)
//       {
//          inta += 64;
//       }
//       if (intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb >= 162 && intb <= 169) 
//       {
//          intb += 64;
//       }
//    } else if (zcode == 218 || zcode == 250) 
//    {
//       /*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code */
//       /*        plus 128 of either lower or upper case 'Z'. */
//       if (inta >= 225 && inta <= 250) 
//       {
//          inta += -32;
//       }
//       if (intb >= 225 && intb <= 250)
//       {
//          intb += -32;
//       }
//    }
//    ret_val = inta == intb;
//    return ret_val;
// } /* lsame_ */

static int
xerbla_ (char *srname, integer * info, ftnlen srname_len)
{
	/* Format strings */
	static char fmt_9999[] = "(\002 ** On entry to \002,a6,\002 parameter nu" "mber \002,i2,\002 had \002,\002an illegal value\002)";

	/* Builtin functions */
	integer s_wsfe (cilist *), do_fio (integer *, char *, ftnlen), e_wsfe (void);
	/* Subroutine */ int s_stop (char *, ftnlen);

	/* Fortran I/O blocks */
	static cilist io___810 = { 0, 6, 0, fmt_9999, 0 };



	/*  -- LAPACK auxiliary routine (version 3.0) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     September 30, 1994 */

	/*     .. Scalar Arguments .. */
	/*     .. */

	/*  Purpose */
	/*  ======= */

	/*  XERBLA  is an error handler for the LAPACK routines. */
	/*  It is called by an LAPACK routine if an input parameter has an */
	/*  invalid value.  A message is printed and execution stops. */

	/*  Installers may consider modifying the STOP statement in order to */
	/*  call system-specific exception-handling facilities. */

	/*  Arguments */
	/*  ========= */

	/*  SRNAME  (input) CHARACTER*6 */
	/*          The name of the routine which called XERBLA. */

	/*  INFO    (input) INTEGER */
	/*          The position of the invalid parameter in the parameter list */
	/*          of the calling routine. */

	/* ===================================================================== */

	/*     .. Executable Statements .. */

	s_wsfe (&io___810);
	do_fio (&c__1, srname, (ftnlen) 6);
	do_fio (&c__1, (char *) &(*info), (ftnlen) sizeof (integer));
	e_wsfe ();

	s_stop ("", (ftnlen) 0);


	/*     End of XERBLA */

	return 0;
}	/* xerbla_ */

/* ----------------------------------------------------------------------- */
/* \BeginDoc */

/* \Name: dgetv0 */

/* \Description: */
/*  Generate a random initial residual vector for the Arnoldi process. */
/*  Force the residual vector to be in the range of the operator OP. */

/* \Usage: */
/*  call dgetv0 */
/*     ( IDO, BMAT, ITRY, INITV, N, J, V, LDV, RESID, RNORM, */
/*       IPNTR, WORKD, IERR ) */

/* \Arguments */
/*  IDO     Integer.  (INPUT/OUTPUT) */
/*          Reverse communication flag.  IDO must be zero on the first */
/*          call to dgetv0. */
/*          ------------------------------------------------------------- */
/*          IDO =  0: first call to the reverse communication interface */
/*          IDO = -1: compute  Y = OP * X  where */
/*                    IPNTR(1) is the pointer into WORKD for X, */
/*                    IPNTR(2) is the pointer into WORKD for Y. */
/*                    This is for the initialization phase to force the */
/*                    starting vector into the range of OP. */
/*          IDO =  2: compute  Y = B * X  where */
/*                    IPNTR(1) is the pointer into WORKD for X, */
/*                    IPNTR(2) is the pointer into WORKD for Y. */
/*          IDO = 99: done */
/*          ------------------------------------------------------------- */

/*  BMAT    Character*1.  (INPUT) */
/*          BMAT specifies the type of the matrix B in the (generalized) */
/*          eigenvalue problem A*x = lambda*B*x. */
/*          B = 'I' -> standard eigenvalue problem A*x = lambda*x */
/*          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x */

/*  ITRY    Integer.  (INPUT) */
/*          ITRY counts the number of times that dgetv0 is called. */
/*          It should be set to 1 on the initial call to dgetv0. */

/*  INITV   Logical variable.  (INPUT) */
/*          .TRUE.  => the initial residual vector is given in RESID. */
/*          .FALSE. => generate a random initial residual vector. */

/*  N       Integer.  (INPUT) */
/*          Dimension of the problem. */

/*  J       Integer.  (INPUT) */
/*          Index of the residual vector to be generated, with respect to */
/*          the Arnoldi process.  J > 1 in case of a "restart". */

/*  V       Double precision N by J array.  (INPUT) */
/*          The first J-1 columns of V contain the current Arnoldi basis */
/*          if this is a "restart". */

/*  LDV     Integer.  (INPUT) */
/*          Leading dimension of V exactly as declared in the calling */
/*          program. */

/*  RESID   Double precision array of length N.  (INPUT/OUTPUT) */
/*          Initial residual vector to be generated.  If RESID is */
/*          provided, force RESID into the range of the operator OP. */

/*  RNORM   Double precision scalar.  (OUTPUT) */
/*          B-norm of the generated residual. */

/*  IPNTR   Integer array of length 3.  (OUTPUT) */

/*  WORKD   Double precision work array of length 2*N.  (REVERSE COMMUNICATION). */
/*          On exit, WORK(1:N) = B*RESID to be used in SSAITR. */

/*  IERR    Integer.  (OUTPUT) */
/*          =  0: Normal exit. */
/*          = -1: Cannot generate a nontrivial restarted residual vector */
/*                in the range of the operator OP. */

/* \EndDoc */

/* ----------------------------------------------------------------------- */

/* \BeginLib */

/* \Local variables: */
/*     xxxxxx  real */

/* \References: */
/*  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in */
/*     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992), */
/*     pp 357-385. */
/*  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly */
/*     Restarted Arnoldi Iteration", Rice University Technical Report */
/*     TR95-13, Department of Computational and Applied Mathematics. */

/* \Routines called: */
/*     second  ARPACK utility routine for timing. */
/*     dvout   ARPACK utility routine for vector output. */
/*     dlarnv  LAPACK routine for generating a random vector. */
/*     dgemv   Level 2 BLAS routine for matrix vector multiplication. */
/*     dcopy   Level 1 BLAS that copies one vector to another. */
/*     ddot    Level 1 BLAS that computes the scalar product of two vectors. */
/*     dnrm2   Level 1 BLAS that computes the norm of a vector. */

/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \SCCS Information: @(#) */
/* FILE: getv0.F   SID: 2.7   DATE OF SID: 04/07/99   RELEASE: 2 */

/* \EndLib */

/* ----------------------------------------------------------------------- */

static int
dgetv0_ (integer * ido, char *bmat, integer * itry, logical
			* initv, integer * n, integer * j, doublereal * v, integer * ldv,
			doublereal * resid, doublereal * rnorm, integer * ipntr, doublereal * workd, integer * ierr, ftnlen bmat_len)
{
	/* Initialized data */

	static logical inits = TRUE_;

	/* System generated locals */
	integer v_dim1, v_offset, i__1;

	/* Builtin functions */
	double sqrt (doublereal);

	/* Local variables */
	static real t0, t1, t2, t3;
	static integer jj;
	static integer iter;
	static logical orth;
	static integer iseed[4];
	static integer idist;
	static logical first;
	static doublereal rnorm0;
	static integer msglvl;


	/*     %----------------------------------------------------% */
	/*     | Include files for debugging and timing information | */
	/*     %----------------------------------------------------% */


	/* \SCCS Information: @(#) */
	/* FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 */

	/*     %---------------------------------% */
	/*     | See debug.doc for documentation | */
	/*     %---------------------------------% */

	/*     %------------------% */
	/*     | Scalar Arguments | */
	/*     %------------------% */

	/*     %--------------------------------% */
	/*     | See stat.doc for documentation | */
	/*     %--------------------------------% */

	/* \SCCS Information: @(#) */
	/* FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 */



	/*     %-----------------% */
	/*     | Array Arguments | */
	/*     %-----------------% */


	/*     %------------% */
	/*     | Parameters | */
	/*     %------------% */


	/*     %------------------------% */
	/*     | Local Scalars & Arrays | */
	/*     %------------------------% */


	/*     %----------------------% */
	/*     | External Subroutines | */
	/*     %----------------------% */


	/*     %--------------------% */
	/*     | External Functions | */
	/*     %--------------------% */


	/*     %---------------------% */
	/*     | Intrinsic Functions | */
	/*     %---------------------% */


	/*     %-----------------% */
	/*     | Data Statements | */
	/*     %-----------------% */

	/* Parameter adjustments */
	--workd;
	--resid;
	v_dim1 = *ldv;
	v_offset = 1 + v_dim1;
	v -= v_offset;
	--ipntr;

	/* Function Body */

	/*     %-----------------------% */
	/*     | Executable Statements | */
	/*     %-----------------------% */


	/*     %-----------------------------------% */
	/*     | Initialize the seed of the LAPACK | */
	/*     | random number generator           | */
	/*     %-----------------------------------% */

	if (inits)
	{
		iseed[0] = 1;
		iseed[1] = 3;
		iseed[2] = 5;
		iseed[3] = 7;
		inits = FALSE_;
	}

	if (*ido == 0)
	{

		/*        %-------------------------------% */
		/*        | Initialize timing statistics  | */
		/*        | & message level for debugging | */
		/*        %-------------------------------% */

		second_ (&t0);
		msglvl = debug_1.mgetv0;

		*ierr = 0;
		iter = 0;
		first = FALSE_;
		orth = FALSE_;

		/*        %-----------------------------------------------------% */
		/*        | Possibly generate a random starting vector in RESID | */
		/*        | Use a LAPACK random number generator used by the    | */
		/*        | matrix generation routines.                         | */
		/*        |    idist = 1: uniform (0,1)  distribution;          | */
		/*        |    idist = 2: uniform (-1,1) distribution;          | */
		/*        |    idist = 3: normal  (0,1)  distribution;          | */
		/*        %-----------------------------------------------------% */

		if (!(*initv))
		{
			idist = 2;
			dlarnv_ (&idist, iseed, n, &resid[1]);
		}

		/*        %----------------------------------------------------------% */
		/*        | Force the starting vector into the range of OP to handle | */
		/*        | the generalized problem when B is possibly (singular).   | */
		/*        %----------------------------------------------------------% */

		second_ (&t2);
		if (*(unsigned char *) bmat == 'G')
		{
			++timing_1.nopx;
			ipntr[1] = 1;
			ipntr[2] = *n + 1;
			dcopy_ (n, &resid[1], &c__1, &workd[1], &c__1);
			*ido = -1;
			goto L9000;
		}
	}

	/*     %-----------------------------------------% */
	/*     | Back from computing OP*(initial-vector) | */
	/*     %-----------------------------------------% */

	if (first)
	{
		goto L20;
	}

	/*     %-----------------------------------------------% */
	/*     | Back from computing B*(orthogonalized-vector) | */
	/*     %-----------------------------------------------% */

	if (orth)
	{
		goto L40;
	}

	if (*(unsigned char *) bmat == 'G')
	{
		second_ (&t3);
		timing_1.tmvopx += t3 - t2;
	}

	/*     %------------------------------------------------------% */
	/*     | Starting vector is now in the range of OP; r = OP*r; | */
	/*     | Compute B-norm of starting vector.                   | */
	/*     %------------------------------------------------------% */

	second_ (&t2);
	first = TRUE_;
	if (*(unsigned char *) bmat == 'G')
	{
		++timing_1.nbx;
		dcopy_ (n, &workd[*n + 1], &c__1, &resid[1], &c__1);
		ipntr[1] = *n + 1;
		ipntr[2] = 1;
		*ido = 2;
		goto L9000;
	}
	else if (*(unsigned char *) bmat == 'I')
	{
		dcopy_ (n, &resid[1], &c__1, &workd[1], &c__1);
	}

 L20:

	if (*(unsigned char *) bmat == 'G')
	{
		second_ (&t3);
		timing_1.tmvbx += t3 - t2;
	}

	first = FALSE_;
	if (*(unsigned char *) bmat == 'G')
	{
		rnorm0 = ddot_ (n, &resid[1], &c__1, &workd[1], &c__1);
		rnorm0 = sqrt ((abs (rnorm0)));
	}
	else if (*(unsigned char *) bmat == 'I')
	{
		rnorm0 = dnrm2_ (n, &resid[1], &c__1);
	}
	*rnorm = rnorm0;

	/*     %---------------------------------------------% */
	/*     | Exit if this is the very first Arnoldi step | */
	/*     %---------------------------------------------% */

	if (*j == 1)
	{
		goto L50;
	}

	/*     %---------------------------------------------------------------- */
	/*     | Otherwise need to B-orthogonalize the starting vector against | */
	/*     | the current Arnoldi basis using Gram-Schmidt with iter. ref.  | */
	/*     | This is the case where an invariant subspace is encountered   | */
	/*     | in the middle of the Arnoldi factorization.                   | */
	/*     |                                                               | */
	/*     |       s = V^{T}*B*r;   r = r - V*s;                           | */
	/*     |                                                               | */
	/*     | Stopping criteria used for iter. ref. is discussed in         | */
	/*     | Parlett's book, page 107 and in Gragg & Reichel TOMS paper.   | */
	/*     %---------------------------------------------------------------% */

	orth = TRUE_;
 L30:

	i__1 = *j - 1;
	dgemv_ ("T", n, &i__1, &c_b348, &v[v_offset], ldv, &workd[1], &c__1, &c_b507, &workd[*n + 1], &c__1, (ftnlen) 1);
	i__1 = *j - 1;
	dgemv_ ("N", n, &i__1, &c_b347, &v[v_offset], ldv, &workd[*n + 1], &c__1, &c_b348, &resid[1], &c__1, (ftnlen) 1);

	/*     %----------------------------------------------------------% */
	/*     | Compute the B-norm of the orthogonalized starting vector | */
	/*     %----------------------------------------------------------% */

	second_ (&t2);
	if (*(unsigned char *) bmat == 'G')
	{
		++timing_1.nbx;
		dcopy_ (n, &resid[1], &c__1, &workd[*n + 1], &c__1);
		ipntr[1] = *n + 1;
		ipntr[2] = 1;
		*ido = 2;
		goto L9000;
	}
	else if (*(unsigned char *) bmat == 'I')
	{
		dcopy_ (n, &resid[1], &c__1, &workd[1], &c__1);
	}

 L40:

	if (*(unsigned char *) bmat == 'G')
	{
		second_ (&t3);
		timing_1.tmvbx += t3 - t2;
	}

	if (*(unsigned char *) bmat == 'G')
	{
		*rnorm = ddot_ (n, &resid[1], &c__1, &workd[1], &c__1);
		*rnorm = sqrt ((abs (*rnorm)));
	}
	else if (*(unsigned char *) bmat == 'I')
	{
		*rnorm = dnrm2_ (n, &resid[1], &c__1);
	}

	/*     %--------------------------------------% */
	/*     | Check for further orthogonalization. | */
	/*     %--------------------------------------% */

	if (msglvl > 2)
	{
		dvout_ (&debug_1.logfil, &c__1, &rnorm0, &debug_1.ndigit, "_getv0: re" "-orthonalization ; rnorm0 is", (ftnlen) 38);
		dvout_ (&debug_1.logfil, &c__1, rnorm, &debug_1.ndigit, "_getv0: re-o" "rthonalization ; rnorm is", (ftnlen) 37);
	}

	if (*rnorm > rnorm0 * .717f)
	{
		goto L50;
	}

	++iter;
	if (iter <= 5)
	{

		/*        %-----------------------------------% */
		/*        | Perform iterative refinement step | */
		/*        %-----------------------------------% */

		rnorm0 = *rnorm;
		goto L30;
	}
	else
	{

		/*        %------------------------------------% */
		/*        | Iterative refinement step "failed" | */
		/*        %------------------------------------% */

		i__1 = *n;
		for (jj = 1; jj <= i__1; ++jj)
		{
			resid[jj] = 0.;
			/* L45: */
		}
		*rnorm = 0.;
		*ierr = -1;
	}

 L50:

	if (msglvl > 0)
	{
		dvout_ (&debug_1.logfil, &c__1, rnorm, &debug_1.ndigit, "_getv0: B-no" "rm of initial / restarted starting vector", (ftnlen) 53);
	}
	if (msglvl > 3)
	{
		dvout_ (&debug_1.logfil, n, &resid[1], &debug_1.ndigit, "_getv0: init" "ial / restarted starting vector", (ftnlen) 43);
	}
	*ido = 99;

	second_ (&t1);
	timing_1.tgetv0 += t1 - t0;

 L9000:
	return 0;

	/*     %---------------% */
	/*     | End of dgetv0 | */
	/*     %---------------% */

}	/* dgetv0_ */

/* ----------------------------------------------------------------------- */
/* \BeginDoc */

/* \Name: dlaqrb */

/* \Description: */
/*  Compute the eigenvalues and the Schur decomposition of an upper */
/*  Hessenberg submatrix in rows and columns ILO to IHI.  Only the */
/*  last component of the Schur vectors are computed. */

/*  This is mostly a modification of the LAPACK routine dlahqr. */

/* \Usage: */
/*  call dlaqrb */
/*     ( WANTT, N, ILO, IHI, H, LDH, WR, WI,  Z, INFO ) */

/* \Arguments */
/*  WANTT   Logical variable.  (INPUT) */
/*          = .TRUE. : the full Schur form T is required; */
/*          = .FALSE.: only eigenvalues are required. */

/*  N       Integer.  (INPUT) */
/*          The order of the matrix H.  N >= 0. */

/*  ILO     Integer.  (INPUT) */
/*  IHI     Integer.  (INPUT) */
/*          It is assumed that H is already upper quasi-triangular in */
/*          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless */
/*          ILO = 1). SLAQRB works primarily with the Hessenberg */
/*          submatrix in rows and columns ILO to IHI, but applies */
/*          transformations to all of H if WANTT is .TRUE.. */
/*          1 <= ILO <= max(1,IHI); IHI <= N. */

/*  H       Double precision array, dimension (LDH,N).  (INPUT/OUTPUT) */
/*          On entry, the upper Hessenberg matrix H. */
/*          On exit, if WANTT is .TRUE., H is upper quasi-triangular in */
/*          rows and columns ILO:IHI, with any 2-by-2 diagonal blocks in */
/*          standard form. If WANTT is .FALSE., the contents of H are */
/*          unspecified on exit. */

/*  LDH     Integer.  (INPUT) */
/*          The leading dimension of the array H. LDH >= max(1,N). */

/*  WR      Double precision array, dimension (N).  (OUTPUT) */
/*  WI      Double precision array, dimension (N).  (OUTPUT) */
/*          The real and imaginary parts, respectively, of the computed */
/*          eigenvalues ILO to IHI are stored in the corresponding */
/*          elements of WR and WI. If two eigenvalues are computed as a */
/*          complex conjugate pair, they are stored in consecutive */
/*          elements of WR and WI, say the i-th and (i+1)th, with */
/*          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the */
/*          eigenvalues are stored in the same order as on the diagonal */
/*          of the Schur form returned in H, with WR(i) = H(i,i), and, if */
/*          H(i:i+1,i:i+1) is a 2-by-2 diagonal block, */
/*          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i). */

/*  Z       Double precision array, dimension (N).  (OUTPUT) */
/*          On exit Z contains the last components of the Schur vectors. */

/*  INFO    Integer.  (OUPUT) */
/*          = 0: successful exit */
/*          > 0: SLAQRB failed to compute all the eigenvalues ILO to IHI */
/*               in a total of 30*(IHI-ILO+1) iterations; if INFO = i, */
/*               elements i+1:ihi of WR and WI contain those eigenvalues */
/*               which have been successfully computed. */

/* \Remarks */
/*  1. None. */

/* ----------------------------------------------------------------------- */

/* \BeginLib */

/* \Local variables: */
/*     xxxxxx  real */

/* \Routines called: */
/*     dlabad  LAPACK routine that computes machine constants. */
/*     dlamch  LAPACK routine that determines machine constants. */
/*     dlanhs  LAPACK routine that computes various norms of a matrix. */
/*     dlanv2  LAPACK routine that computes the Schur factorization of */
/*             2 by 2 nonsymmetric matrix in standard form. */
/*     dlarfg  LAPACK Householder reflection construction routine. */
/*     dcopy   Level 1 BLAS that copies one vector to another. */
/*     drot    Level 1 BLAS that applies a rotation to a 2 by 2 matrix. */

/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \Revision history: */
/*     xx/xx/92: Version ' 2.4' */
/*               Modified from the LAPACK routine dlahqr so that only the */
/*               last component of the Schur vectors are computed. */

/* \SCCS Information: @(#) */
/* FILE: laqrb.F   SID: 2.2   DATE OF SID: 8/27/96   RELEASE: 2 */

/* \Remarks */
/*     1. None */

/* \EndLib */

/* ----------------------------------------------------------------------- */

static int
dlaqrb_ (logical * wantt, integer * n, integer * ilo,
			integer * ihi, doublereal * h__, integer * ldh, doublereal * wr, doublereal * wi, doublereal * z__, integer * info)
{
	/* System generated locals */
	integer h_dim1, h_offset, i__1, i__2, i__3, i__4;
	doublereal d__1, d__2;

	/* Local variables */
	static integer i__, j, k, l, m;
	static doublereal s, v[3];
	static integer i1, i2;
	static doublereal t1, t2, t3, v1, v2, v3, h00, h10, h11, h12, h21, h22, h33, h44;
	static integer nh;
	static doublereal cs;
	static integer nr;
	static doublereal sn, h33s, h44s;
	static integer itn, its;
	static doublereal ulp, sum, tst1, h43h34, unfl, ovfl;
	static doublereal work[1];
	static doublereal smlnum;


	/*     %------------------% */
	/*     | Scalar Arguments | */
	/*     %------------------% */


	/*     %-----------------% */
	/*     | Array Arguments | */
	/*     %-----------------% */


	/*     %------------% */
	/*     | Parameters | */
	/*     %------------% */


	/*     %------------------------% */
	/*     | Local Scalars & Arrays | */
	/*     %------------------------% */


	/*     %--------------------% */
	/*     | External Functions | */
	/*     %--------------------% */


	/*     %----------------------% */
	/*     | External Subroutines | */
	/*     %----------------------% */


	/*     %-----------------------% */
	/*     | Executable Statements | */
	/*     %-----------------------% */

	/* Parameter adjustments */
	h_dim1 = *ldh;
	h_offset = 1 + h_dim1;
	h__ -= h_offset;
	--wr;
	--wi;
	--z__;

	/* Function Body */
	*info = 0;

	/*     %--------------------------% */
	/*     | Quick return if possible | */
	/*     %--------------------------% */

	if (*n == 0)
	{
		return 0;
	}
	if (*ilo == *ihi)
	{
		wr[*ilo] = h__[*ilo + *ilo * h_dim1];
		wi[*ilo] = 0.;
		return 0;
	}

	/*     %---------------------------------------------% */
	/*     | Initialize the vector of last components of | */
	/*     | the Schur vectors for accumulation.         | */
	/*     %---------------------------------------------% */

	i__1 = *n - 1;
	for (j = 1; j <= i__1; ++j)
	{
		z__[j] = 0.;
		/* L5: */
	}
	z__[*n] = 1.;

	nh = *ihi - *ilo + 1;

	/*     %-------------------------------------------------------------% */
	/*     | Set machine-dependent constants for the stopping criterion. | */
	/*     | If norm(H) <= sqrt(OVFL), overflow should not occur.        | */
	/*     %-------------------------------------------------------------% */

	unfl = dlamch_ ("safe minimum", (ftnlen) 12);
	ovfl = 1. / unfl;
	dlabad_ (&unfl, &ovfl);
	ulp = dlamch_ ("precision", (ftnlen) 9);
	smlnum = unfl * (nh / ulp);

	/*     %---------------------------------------------------------------% */
	/*     | I1 and I2 are the indices of the first row and last column    | */
	/*     | of H to which transformations must be applied. If eigenvalues | */
	/*     | only are computed, I1 and I2 are set inside the main loop.    | */
	/*     | Zero out H(J+2,J) = ZERO for J=1:N if WANTT = .TRUE.          | */
	/*     | else H(J+2,J) for J=ILO:IHI-ILO-1 if WANTT = .FALSE.          | */
	/*     %---------------------------------------------------------------% */

	if (*wantt)
	{
		i1 = 1;
		i2 = *n;
		i__1 = i2 - 2;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			h__[i1 + i__ + 1 + i__ * h_dim1] = 0.;
			/* L8: */
		}
	}
	else
	{
		i__1 = *ihi - *ilo - 1;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			h__[*ilo + i__ + 1 + (*ilo + i__ - 1) * h_dim1] = 0.;
			/* L9: */
		}
	}

	/*     %---------------------------------------------------% */
	/*     | ITN is the total number of QR iterations allowed. | */
	/*     %---------------------------------------------------% */

	itn = nh * 30;

	/*     ------------------------------------------------------------------ */
	/*     The main loop begins here. I is the loop index and decreases from */
	/*     IHI to ILO in steps of 1 or 2. Each iteration of the loop works */
	/*     with the active submatrix in rows and columns L to I. */
	/*     Eigenvalues I+1 to IHI have already converged. Either L = ILO or */
	/*     H(L,L-1) is negligible so that the matrix splits. */
	/*     ------------------------------------------------------------------ */

	i__ = *ihi;
 L10:
	l = *ilo;
	if (i__ < *ilo)
	{
		goto L150;
	}
	/*     %--------------------------------------------------------------% */
	/*     | Perform QR iterations on rows and columns ILO to I until a   | */
	/*     | submatrix of order 1 or 2 splits off at the bottom because a | */
	/*     | subdiagonal element has become negligible.                   | */
	/*     %--------------------------------------------------------------% */
	i__1 = itn;
	for (its = 0; its <= i__1; ++its)
	{

		/*        %----------------------------------------------% */
		/*        | Look for a single small subdiagonal element. | */
		/*        %----------------------------------------------% */

		i__2 = l + 1;
		for (k = i__; k >= i__2; --k)
		{
			tst1 = (d__1 = h__[k - 1 + (k - 1) * h_dim1], abs (d__1)) + (d__2 = h__[k + k * h_dim1], abs (d__2));
			if (tst1 == 0.)
			{
				i__3 = i__ - l + 1;
				tst1 = dlanhs_ ("1", &i__3, &h__[l + l * h_dim1], ldh, work, (ftnlen) 1);
			}
			/* Computing MAX */
			d__2 = ulp * tst1;
			if ((d__1 = h__[k + (k - 1) * h_dim1], abs (d__1)) <= max (d__2, smlnum))
			{
				goto L30;
			}
			/* L20: */
		}
	 L30:
		l = k;
		if (l > *ilo)
		{

			/*           %------------------------% */
			/*           | H(L,L-1) is negligible | */
			/*           %------------------------% */

			h__[l + (l - 1) * h_dim1] = 0.;
		}

		/*        %-------------------------------------------------------------% */
		/*        | Exit from loop if a submatrix of order 1 or 2 has split off | */
		/*        %-------------------------------------------------------------% */

		if (l >= i__ - 1)
		{
			goto L140;
		}

		/*        %---------------------------------------------------------% */
		/*        | Now the active submatrix is in rows and columns L to I. | */
		/*        | If eigenvalues only are being computed, only the active | */
		/*        | submatrix need be transformed.                          | */
		/*        %---------------------------------------------------------% */

		if (!(*wantt))
		{
			i1 = l;
			i2 = i__;
		}

		if (its == 10 || its == 20)
		{

			/*           %-------------------% */
			/*           | Exceptional shift | */
			/*           %-------------------% */

			s = (d__1 = h__[i__ + (i__ - 1) * h_dim1], abs (d__1)) + (d__2 = h__[i__ - 1 + (i__ - 2) * h_dim1], abs (d__2));
			h44 = s * .75;
			h33 = h44;
			h43h34 = s * -.4375 * s;

		}
		else
		{

			/*           %-----------------------------------------% */
			/*           | Prepare to use Wilkinson's double shift | */
			/*           %-----------------------------------------% */

			h44 = h__[i__ + i__ * h_dim1];
			h33 = h__[i__ - 1 + (i__ - 1) * h_dim1];
			h43h34 = h__[i__ + (i__ - 1) * h_dim1] * h__[i__ - 1 + i__ * h_dim1];
		}

		/*        %-----------------------------------------------------% */
		/*        | Look for two consecutive small subdiagonal elements | */
		/*        %-----------------------------------------------------% */

		i__2 = l;
		for (m = i__ - 2; m >= i__2; --m)
		{

			/*           %---------------------------------------------------------% */
			/*           | Determine the effect of starting the double-shift QR    | */
			/*           | iteration at row M, and see if this would make H(M,M-1) | */
			/*           | negligible.                                             | */
			/*           %---------------------------------------------------------% */

			h11 = h__[m + m * h_dim1];
			h22 = h__[m + 1 + (m + 1) * h_dim1];
			h21 = h__[m + 1 + m * h_dim1];
			h12 = h__[m + (m + 1) * h_dim1];
			h44s = h44 - h11;
			h33s = h33 - h11;
			v1 = (h33s * h44s - h43h34) / h21 + h12;
			v2 = h22 - h11 - h33s - h44s;
			v3 = h__[m + 2 + (m + 1) * h_dim1];
			s = abs (v1) + abs (v2) + abs (v3);
			v1 /= s;
			v2 /= s;
			v3 /= s;
			v[0] = v1;
			v[1] = v2;
			v[2] = v3;
			if (m == l)
			{
				goto L50;
			}
			h00 = h__[m - 1 + (m - 1) * h_dim1];
			h10 = h__[m + (m - 1) * h_dim1];
			tst1 = abs (v1) * (abs (h00) + abs (h11) + abs (h22));
			if (abs (h10) * (abs (v2) + abs (v3)) <= ulp * tst1)
			{
				goto L50;
			}
			/* L40: */
		}
	 L50:

		/*        %----------------------% */
		/*        | Double-shift QR step | */
		/*        %----------------------% */

		i__2 = i__ - 1;
		for (k = m; k <= i__2; ++k)
		{

			/*           ------------------------------------------------------------ */
			/*           The first iteration of this loop determines a reflection G */
			/*           from the vector V and applies it from left and right to H, */
			/*           thus creating a nonzero bulge below the subdiagonal. */

			/*           Each subsequent iteration determines a reflection G to */
			/*           restore the Hessenberg form in the (K-1)th column, and thus */
			/*           chases the bulge one step toward the bottom of the active */
			/*           submatrix. NR is the order of G. */
			/*           ------------------------------------------------------------ */

			/* Computing MIN */
			i__3 = 3, i__4 = i__ - k + 1;
			nr = min (i__3, i__4);
			if (k > m)
			{
				dcopy_ (&nr, &h__[k + (k - 1) * h_dim1], &c__1, v, &c__1);
			}
			dlarfg_ (&nr, v, &v[1], &c__1, &t1);
			if (k > m)
			{
				h__[k + (k - 1) * h_dim1] = v[0];
				h__[k + 1 + (k - 1) * h_dim1] = 0.;
				if (k < i__ - 1)
				{
					h__[k + 2 + (k - 1) * h_dim1] = 0.;
				}
			}
			else if (m > l)
			{
				h__[k + (k - 1) * h_dim1] = -h__[k + (k - 1) * h_dim1];
			}
			v2 = v[1];
			t2 = t1 * v2;
			if (nr == 3)
			{
				v3 = v[2];
				t3 = t1 * v3;

				/*              %------------------------------------------------% */
				/*              | Apply G from the left to transform the rows of | */
				/*              | the matrix in columns K to I2.                 | */
				/*              %------------------------------------------------% */

				i__3 = i2;
				for (j = k; j <= i__3; ++j)
				{
					sum = h__[k + j * h_dim1] + v2 * h__[k + 1 + j * h_dim1] + v3 * h__[k + 2 + j * h_dim1];
					h__[k + j * h_dim1] -= sum * t1;
					h__[k + 1 + j * h_dim1] -= sum * t2;
					h__[k + 2 + j * h_dim1] -= sum * t3;
					/* L60: */
				}

				/*              %----------------------------------------------------% */
				/*              | Apply G from the right to transform the columns of | */
				/*              | the matrix in rows I1 to min(K+3,I).               | */
				/*              %----------------------------------------------------% */

				/* Computing MIN */
				i__4 = k + 3;
				i__3 = min (i__4, i__);
				for (j = i1; j <= i__3; ++j)
				{
					sum = h__[j + k * h_dim1] + v2 * h__[j + (k + 1) * h_dim1] + v3 * h__[j + (k + 2) * h_dim1];
					h__[j + k * h_dim1] -= sum * t1;
					h__[j + (k + 1) * h_dim1] -= sum * t2;
					h__[j + (k + 2) * h_dim1] -= sum * t3;
					/* L70: */
				}

				/*              %----------------------------------% */
				/*              | Accumulate transformations for Z | */
				/*              %----------------------------------% */

				sum = z__[k] + v2 * z__[k + 1] + v3 * z__[k + 2];
				z__[k] -= sum * t1;
				z__[k + 1] -= sum * t2;
				z__[k + 2] -= sum * t3;
			}
			else if (nr == 2)
			{

				/*              %------------------------------------------------% */
				/*              | Apply G from the left to transform the rows of | */
				/*              | the matrix in columns K to I2.                 | */
				/*              %------------------------------------------------% */

				i__3 = i2;
				for (j = k; j <= i__3; ++j)
				{
					sum = h__[k + j * h_dim1] + v2 * h__[k + 1 + j * h_dim1];
					h__[k + j * h_dim1] -= sum * t1;
					h__[k + 1 + j * h_dim1] -= sum * t2;
					/* L90: */
				}

				/*              %----------------------------------------------------% */
				/*              | Apply G from the right to transform the columns of | */
				/*              | the matrix in rows I1 to min(K+3,I).               | */
				/*              %----------------------------------------------------% */

				i__3 = i__;
				for (j = i1; j <= i__3; ++j)
				{
					sum = h__[j + k * h_dim1] + v2 * h__[j + (k + 1) * h_dim1];
					h__[j + k * h_dim1] -= sum * t1;
					h__[j + (k + 1) * h_dim1] -= sum * t2;
					/* L100: */
				}

				/*              %----------------------------------% */
				/*              | Accumulate transformations for Z | */
				/*              %----------------------------------% */

				sum = z__[k] + v2 * z__[k + 1];
				z__[k] -= sum * t1;
				z__[k + 1] -= sum * t2;
			}
			/* L120: */
		}
		/* L130: */
	}

	/*     %-------------------------------------------------------% */
	/*     | Failure to converge in remaining number of iterations | */
	/*     %-------------------------------------------------------% */

	*info = i__;
	return 0;
 L140:
	if (l == i__)
	{

		/*        %------------------------------------------------------% */
		/*        | H(I,I-1) is negligible: one eigenvalue has converged | */
		/*        %------------------------------------------------------% */

		wr[i__] = h__[i__ + i__ * h_dim1];
		wi[i__] = 0.;
	}
	else if (l == i__ - 1)
	{

		/*        %--------------------------------------------------------% */
		/*        | H(I-1,I-2) is negligible;                              | */
		/*        | a pair of eigenvalues have converged.                  | */
		/*        |                                                        | */
		/*        | Transform the 2-by-2 submatrix to standard Schur form, | */
		/*        | and compute and store the eigenvalues.                 | */
		/*        %--------------------------------------------------------% */

		dlanv2_ (&h__[i__ - 1 + (i__ - 1) * h_dim1], &h__[i__ - 1 + i__ *
																		  h_dim1], &h__[i__ + (i__ - 1) * h_dim1], &h__[i__ + i__ *
																																		h_dim1], &wr[i__ - 1], &wi[i__ - 1], &wr[i__], &wi[i__],
					&cs, &sn);
		if (*wantt)
		{

			/*           %-----------------------------------------------------% */
			/*           | Apply the transformation to the rest of H and to Z, | */
			/*           | as required.                                        | */
			/*           %-----------------------------------------------------% */

			if (i2 > i__)
			{
				i__1 = i2 - i__;
				drot_ (&i__1, &h__[i__ - 1 + (i__ + 1) * h_dim1], ldh, &h__[i__ + (i__ + 1) * h_dim1], ldh, &cs, &sn);
			}
			i__1 = i__ - i1 - 1;
			drot_ (&i__1, &h__[i1 + (i__ - 1) * h_dim1], &c__1, &h__[i1 + i__ * h_dim1], &c__1, &cs, &sn);
			sum = cs * z__[i__ - 1] + sn * z__[i__];
			z__[i__] = cs * z__[i__] - sn * z__[i__ - 1];
			z__[i__ - 1] = sum;
		}
	}

	/*     %---------------------------------------------------------% */
	/*     | Decrement number of remaining iterations, and return to | */
	/*     | start of the main loop with new value of I.             | */
	/*     %---------------------------------------------------------% */

	itn -= its;
	i__ = l - 1;
	goto L10;
 L150:
	return 0;

	/*     %---------------% */
	/*     | End of dlaqrb | */
	/*     %---------------% */

}	/* dlaqrb_ */

/* ----------------------------------------------------------------------- */
/*  Routine:    DMOUT */

/*  Purpose:    Real matrix output routine. */

/*  Usage:      CALL DMOUT (LOUT, M, N, A, LDA, IDIGIT, IFMT) */

/*  Arguments */
/*     M      - Number of rows of A.  (Input) */
/*     N      - Number of columns of A.  (Input) */
/*     A      - Real M by N matrix to be printed.  (Input) */
/*     LDA    - Leading dimension of A exactly as specified in the */
/*              dimension statement of the calling program.  (Input) */
/*     IFMT   - Format to be used in printing matrix A.  (Input) */
/*     IDIGIT - Print up to IABS(IDIGIT) decimal digits per number.  (In) */
/*              If IDIGIT .LT. 0, printing is done with 72 columns. */
/*              If IDIGIT .GT. 0, printing is done with 132 columns. */

/* ----------------------------------------------------------------------- */

static int
dmout_ (integer * lout, integer * m, integer * n, doublereal * a, integer * lda, integer * idigit, char *ifmt, ftnlen ifmt_len)
{
	/* Initialized data */

	static char icol[1 * 3] = "C" "o" "l";

	/* Format strings */
	static char fmt_9999[] = "(/1x,a,/1x,a)";
	static char fmt_9998[] = "(10x,10(4x,3a1,i4,1x))";
	static char fmt_9994[] = "(1x,\002 Row\002,i4,\002:\002,1x,1p,10d12.3)";
	static char fmt_9997[] = "(10x,8(5x,3a1,i4,2x))";
	static char fmt_9993[] = "(1x,\002 Row\002,i4,\002:\002,1x,1p,8d14.5)";
	static char fmt_9996[] = "(10x,6(7x,3a1,i4,4x))";
	static char fmt_9992[] = "(1x,\002 Row\002,i4,\002:\002,1x,1p,6d18.9)";
	static char fmt_9995[] = "(10x,5(9x,3a1,i4,6x))";
	static char fmt_9991[] = "(1x,\002 Row\002,i4,\002:\002,1x,1p,5d22.13)";
	static char fmt_9990[] = "(1x,\002 \002)";

	/* System generated locals */
	integer a_dim1, a_offset, i__1, i__2, i__3;

	/* Builtin functions */
	integer i_len (char *, ftnlen), s_wsfe (cilist *), do_fio (integer *, char *, ftnlen), e_wsfe (void);

	/* Local variables */
	static integer i__, j, k1, k2, lll;
	static char line[80];
	static integer ndigit;

	/* Fortran I/O blocks */
	static cilist io___867 = { 0, 0, 0, fmt_9999, 0 };
	static cilist io___871 = { 0, 0, 0, fmt_9998, 0 };
	static cilist io___872 = { 0, 0, 0, fmt_9994, 0 };
	static cilist io___874 = { 0, 0, 0, fmt_9997, 0 };
	static cilist io___875 = { 0, 0, 0, fmt_9993, 0 };
	static cilist io___876 = { 0, 0, 0, fmt_9996, 0 };
	static cilist io___877 = { 0, 0, 0, fmt_9992, 0 };
	static cilist io___878 = { 0, 0, 0, fmt_9995, 0 };
	static cilist io___879 = { 0, 0, 0, fmt_9991, 0 };
	static cilist io___880 = { 0, 0, 0, fmt_9998, 0 };
	static cilist io___881 = { 0, 0, 0, fmt_9994, 0 };
	static cilist io___882 = { 0, 0, 0, fmt_9997, 0 };
	static cilist io___883 = { 0, 0, 0, fmt_9993, 0 };
	static cilist io___884 = { 0, 0, 0, fmt_9996, 0 };
	static cilist io___885 = { 0, 0, 0, fmt_9992, 0 };
	static cilist io___886 = { 0, 0, 0, fmt_9995, 0 };
	static cilist io___887 = { 0, 0, 0, fmt_9991, 0 };
	static cilist io___888 = { 0, 0, 0, fmt_9990, 0 };


	/*     ... */
	/*     ... SPECIFICATIONS FOR ARGUMENTS */
	/*     ... */
	/*     ... SPECIFICATIONS FOR LOCAL VARIABLES */
	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Local Arrays .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Data statements .. */
	/* Parameter adjustments */
	a_dim1 = *lda;
	a_offset = 1 + a_dim1;
	a -= a_offset;

	/* Function Body */
	/*     .. */
	/*     .. Executable Statements .. */
	/*     ... */
	/*     ... FIRST EXECUTABLE STATEMENT */

	/* Computing MIN */
	i__1 = i_len (ifmt, ifmt_len);
	lll = min (i__1, 80);
	i__1 = lll;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		*(unsigned char *) &line[i__ - 1] = '-';
		/* L10: */
	}

	for (i__ = lll + 1; i__ <= 80; ++i__)
	{
		*(unsigned char *) &line[i__ - 1] = ' ';
		/* L20: */
	}

	io___867.ciunit = *lout;
	s_wsfe (&io___867);
	do_fio (&c__1, ifmt, ifmt_len);
	do_fio (&c__1, line, lll);
	e_wsfe ();

	if (*m <= 0 || *n <= 0 || *lda <= 0)
	{
		return 0;
	}
	ndigit = *idigit;
	if (*idigit == 0)
	{
		ndigit = 4;
	}

	/* ======================================================================= */
	/*             CODE FOR OUTPUT USING 72 COLUMNS FORMAT */
	/* ======================================================================= */

	if (*idigit < 0)
	{
		ndigit = -(*idigit);
		if (ndigit <= 4)
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 5)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 4;
				k2 = min (i__2, i__3);
				io___871.ciunit = *lout;
				s_wsfe (&io___871);
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__3, icol, (ftnlen) 1);
					do_fio (&c__1, (char *) &i__, (ftnlen) sizeof (integer));
				}
				e_wsfe ();
				i__2 = *m;
				for (i__ = 1; i__ <= i__2; ++i__)
				{
					io___872.ciunit = *lout;
					s_wsfe (&io___872);
					do_fio (&c__1, (char *) &i__, (ftnlen) sizeof (integer));
					i__3 = k2;
					for (j = k1; j <= i__3; ++j)
					{
						do_fio (&c__1, (char *) &a[i__ + j * a_dim1], (ftnlen) sizeof (doublereal));
					}
					e_wsfe ();
					/* L30: */
				}
				/* L40: */
			}

		}
		else if (ndigit <= 6)
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 4)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 3;
				k2 = min (i__2, i__3);
				io___874.ciunit = *lout;
				s_wsfe (&io___874);
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__3, icol, (ftnlen) 1);
					do_fio (&c__1, (char *) &i__, (ftnlen) sizeof (integer));
				}
				e_wsfe ();
				i__2 = *m;
				for (i__ = 1; i__ <= i__2; ++i__)
				{
					io___875.ciunit = *lout;
					s_wsfe (&io___875);
					do_fio (&c__1, (char *) &i__, (ftnlen) sizeof (integer));
					i__3 = k2;
					for (j = k1; j <= i__3; ++j)
					{
						do_fio (&c__1, (char *) &a[i__ + j * a_dim1], (ftnlen) sizeof (doublereal));
					}
					e_wsfe ();
					/* L50: */
				}
				/* L60: */
			}

		}
		else if (ndigit <= 10)
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 3)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 2;
				k2 = min (i__2, i__3);
				io___876.ciunit = *lout;
				s_wsfe (&io___876);
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__3, icol, (ftnlen) 1);
					do_fio (&c__1, (char *) &i__, (ftnlen) sizeof (integer));
				}
				e_wsfe ();
				i__2 = *m;
				for (i__ = 1; i__ <= i__2; ++i__)
				{
					io___877.ciunit = *lout;
					s_wsfe (&io___877);
					do_fio (&c__1, (char *) &i__, (ftnlen) sizeof (integer));
					i__3 = k2;
					for (j = k1; j <= i__3; ++j)
					{
						do_fio (&c__1, (char *) &a[i__ + j * a_dim1], (ftnlen) sizeof (doublereal));
					}
					e_wsfe ();
					/* L70: */
				}
				/* L80: */
			}

		}
		else
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 2)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 1;
				k2 = min (i__2, i__3);
				io___878.ciunit = *lout;
				s_wsfe (&io___878);
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__3, icol, (ftnlen) 1);
					do_fio (&c__1, (char *) &i__, (ftnlen) sizeof (integer));
				}
				e_wsfe ();
				i__2 = *m;
				for (i__ = 1; i__ <= i__2; ++i__)
				{
					io___879.ciunit = *lout;
					s_wsfe (&io___879);
					do_fio (&c__1, (char *) &i__, (ftnlen) sizeof (integer));
					i__3 = k2;
					for (j = k1; j <= i__3; ++j)
					{
						do_fio (&c__1, (char *) &a[i__ + j * a_dim1], (ftnlen) sizeof (doublereal));
					}
					e_wsfe ();
					/* L90: */
				}
				/* L100: */
			}
		}

		/* ======================================================================= */
		/*             CODE FOR OUTPUT USING 132 COLUMNS FORMAT */
		/* ======================================================================= */

	}
	else
	{
		if (ndigit <= 4)
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 10)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 9;
				k2 = min (i__2, i__3);
				io___880.ciunit = *lout;
				s_wsfe (&io___880);
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__3, icol, (ftnlen) 1);
					do_fio (&c__1, (char *) &i__, (ftnlen) sizeof (integer));
				}
				e_wsfe ();
				i__2 = *m;
				for (i__ = 1; i__ <= i__2; ++i__)
				{
					io___881.ciunit = *lout;
					s_wsfe (&io___881);
					do_fio (&c__1, (char *) &i__, (ftnlen) sizeof (integer));
					i__3 = k2;
					for (j = k1; j <= i__3; ++j)
					{
						do_fio (&c__1, (char *) &a[i__ + j * a_dim1], (ftnlen) sizeof (doublereal));
					}
					e_wsfe ();
					/* L110: */
				}
				/* L120: */
			}

		}
		else if (ndigit <= 6)
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 8)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 7;
				k2 = min (i__2, i__3);
				io___882.ciunit = *lout;
				s_wsfe (&io___882);
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__3, icol, (ftnlen) 1);
					do_fio (&c__1, (char *) &i__, (ftnlen) sizeof (integer));
				}
				e_wsfe ();
				i__2 = *m;
				for (i__ = 1; i__ <= i__2; ++i__)
				{
					io___883.ciunit = *lout;
					s_wsfe (&io___883);
					do_fio (&c__1, (char *) &i__, (ftnlen) sizeof (integer));
					i__3 = k2;
					for (j = k1; j <= i__3; ++j)
					{
						do_fio (&c__1, (char *) &a[i__ + j * a_dim1], (ftnlen) sizeof (doublereal));
					}
					e_wsfe ();
					/* L130: */
				}
				/* L140: */
			}

		}
		else if (ndigit <= 10)
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 6)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 5;
				k2 = min (i__2, i__3);
				io___884.ciunit = *lout;
				s_wsfe (&io___884);
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__3, icol, (ftnlen) 1);
					do_fio (&c__1, (char *) &i__, (ftnlen) sizeof (integer));
				}
				e_wsfe ();
				i__2 = *m;
				for (i__ = 1; i__ <= i__2; ++i__)
				{
					io___885.ciunit = *lout;
					s_wsfe (&io___885);
					do_fio (&c__1, (char *) &i__, (ftnlen) sizeof (integer));
					i__3 = k2;
					for (j = k1; j <= i__3; ++j)
					{
						do_fio (&c__1, (char *) &a[i__ + j * a_dim1], (ftnlen) sizeof (doublereal));
					}
					e_wsfe ();
					/* L150: */
				}
				/* L160: */
			}

		}
		else
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 5)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 4;
				k2 = min (i__2, i__3);
				io___886.ciunit = *lout;
				s_wsfe (&io___886);
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__3, icol, (ftnlen) 1);
					do_fio (&c__1, (char *) &i__, (ftnlen) sizeof (integer));
				}
				e_wsfe ();
				i__2 = *m;
				for (i__ = 1; i__ <= i__2; ++i__)
				{
					io___887.ciunit = *lout;
					s_wsfe (&io___887);
					do_fio (&c__1, (char *) &i__, (ftnlen) sizeof (integer));
					i__3 = k2;
					for (j = k1; j <= i__3; ++j)
					{
						do_fio (&c__1, (char *) &a[i__ + j * a_dim1], (ftnlen) sizeof (doublereal));
					}
					e_wsfe ();
					/* L170: */
				}
				/* L180: */
			}
		}
	}
	io___888.ciunit = *lout;
	s_wsfe (&io___888);
	e_wsfe ();


	return 0;
}	/* dmout_ */

/* ----------------------------------------------------------------------- */
/* \BeginDoc */

/* \Name: dnaitr */

/* \Description: */
/*  Reverse communication interface for applying NP additional steps to */
/*  a K step nonsymmetric Arnoldi factorization. */

/*  Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T */

/*          with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0. */

/*  Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T */

/*          with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0. */

/*  where OP and B are as in dnaupd.  The B-norm of r_{k+p} is also */
/*  computed and returned. */

/* \Usage: */
/*  call dnaitr */
/*     ( IDO, BMAT, N, K, NP, NB, RESID, RNORM, V, LDV, H, LDH, */
/*       IPNTR, WORKD, INFO ) */

/* \Arguments */
/*  IDO     Integer.  (INPUT/OUTPUT) */
/*          Reverse communication flag. */
/*          ------------------------------------------------------------- */
/*          IDO =  0: first call to the reverse communication interface */
/*          IDO = -1: compute  Y = OP * X  where */
/*                    IPNTR(1) is the pointer into WORK for X, */
/*                    IPNTR(2) is the pointer into WORK for Y. */
/*                    This is for the restart phase to force the new */
/*                    starting vector into the range of OP. */
/*          IDO =  1: compute  Y = OP * X  where */
/*                    IPNTR(1) is the pointer into WORK for X, */
/*                    IPNTR(2) is the pointer into WORK for Y, */
/*                    IPNTR(3) is the pointer into WORK for B * X. */
/*          IDO =  2: compute  Y = B * X  where */
/*                    IPNTR(1) is the pointer into WORK for X, */
/*                    IPNTR(2) is the pointer into WORK for Y. */
/*          IDO = 99: done */
/*          ------------------------------------------------------------- */
/*          When the routine is used in the "shift-and-invert" mode, the */
/*          vector B * Q is already available and do not need to be */
/*          recompute in forming OP * Q. */

/*  BMAT    Character*1.  (INPUT) */
/*          BMAT specifies the type of the matrix B that defines the */
/*          semi-inner product for the operator OP.  See dnaupd. */
/*          B = 'I' -> standard eigenvalue problem A*x = lambda*x */
/*          B = 'G' -> generalized eigenvalue problem A*x = lambda*M**x */

/*  N       Integer.  (INPUT) */
/*          Dimension of the eigenproblem. */

/*  K       Integer.  (INPUT) */
/*          Current size of V and H. */

/*  NP      Integer.  (INPUT) */
/*          Number of additional Arnoldi steps to take. */

/*  NB      Integer.  (INPUT) */
/*          Blocksize to be used in the recurrence. */
/*          Only work for NB = 1 right now.  The goal is to have a */
/*          program that implement both the block and non-block method. */

/*  RESID   Double precision array of length N.  (INPUT/OUTPUT) */
/*          On INPUT:  RESID contains the residual vector r_{k}. */
/*          On OUTPUT: RESID contains the residual vector r_{k+p}. */

/*  RNORM   Double precision scalar.  (INPUT/OUTPUT) */
/*          B-norm of the starting residual on input. */
/*          B-norm of the updated residual r_{k+p} on output. */

/*  V       Double precision N by K+NP array.  (INPUT/OUTPUT) */
/*          On INPUT:  V contains the Arnoldi vectors in the first K */
/*          columns. */
/*          On OUTPUT: V contains the new NP Arnoldi vectors in the next */
/*          NP columns.  The first K columns are unchanged. */

/*  LDV     Integer.  (INPUT) */
/*          Leading dimension of V exactly as declared in the calling */
/*          program. */

/*  H       Double precision (K+NP) by (K+NP) array.  (INPUT/OUTPUT) */
/*          H is used to store the generated upper Hessenberg matrix. */

/*  LDH     Integer.  (INPUT) */
/*          Leading dimension of H exactly as declared in the calling */
/*          program. */

/*  IPNTR   Integer array of length 3.  (OUTPUT) */
/*          Pointer to mark the starting locations in the WORK for */
/*          vectors used by the Arnoldi iteration. */
/*          ------------------------------------------------------------- */
/*          IPNTR(1): pointer to the current operand vector X. */
/*          IPNTR(2): pointer to the current result vector Y. */
/*          IPNTR(3): pointer to the vector B * X when used in the */
/*                    shift-and-invert mode.  X is the current operand. */
/*          ------------------------------------------------------------- */

/*  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION) */
/*          Distributed array to be used in the basic Arnoldi iteration */
/*          for reverse communication.  The calling program should not */
/*          use WORKD as temporary workspace during the iteration !!!!!! */
/*          On input, WORKD(1:N) = B*RESID and is used to save some */
/*          computation at the first step. */

/*  INFO    Integer.  (OUTPUT) */
/*          = 0: Normal exit. */
/*          > 0: Size of the spanning invariant subspace of OP found. */

/* \EndDoc */

/* ----------------------------------------------------------------------- */

/* \BeginLib */

/* \Local variables: */
/*     xxxxxx  real */

/* \References: */
/*  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in */
/*     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992), */
/*     pp 357-385. */
/*  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly */
/*     Restarted Arnoldi Iteration", Rice University Technical Report */
/*     TR95-13, Department of Computational and Applied Mathematics. */

/* \Routines called: */
/*     dgetv0  ARPACK routine to generate the initial vector. */
/*     ivout   ARPACK utility routine that prints integers. */
/*     second  ARPACK utility routine for timing. */
/*     dmout   ARPACK utility routine that prints matrices */
/*     dvout   ARPACK utility routine that prints vectors. */
/*     dlabad  LAPACK routine that computes machine constants. */
/*     dlamch  LAPACK routine that determines machine constants. */
/*     dlascl  LAPACK routine for careful scaling of a matrix. */
/*     dlanhs  LAPACK routine that computes various norms of a matrix. */
/*     dgemv   Level 2 BLAS routine for matrix vector multiplication. */
/*     daxpy   Level 1 BLAS that computes a vector triad. */
/*     dscal   Level 1 BLAS that scales a vector. */
/*     dcopy   Level 1 BLAS that copies one vector to another . */
/*     ddot    Level 1 BLAS that computes the scalar product of two vectors. */
/*     dnrm2   Level 1 BLAS that computes the norm of a vector. */

/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \Revision history: */
/*     xx/xx/92: Version ' 2.4' */

/* \SCCS Information: @(#) */
/* FILE: naitr.F   SID: 2.4   DATE OF SID: 8/27/96   RELEASE: 2 */

/* \Remarks */
/*  The algorithm implemented is: */

/*  restart = .false. */
/*  Given V_{k} = [v_{1}, ..., v_{k}], r_{k}; */
/*  r_{k} contains the initial residual vector even for k = 0; */
/*  Also assume that rnorm = || B*r_{k} || and B*r_{k} are already */
/*  computed by the calling program. */

/*  betaj = rnorm ; p_{k+1} = B*r_{k} ; */
/*  For  j = k+1, ..., k+np  Do */
/*     1) if ( betaj < tol ) stop or restart depending on j. */
/*        ( At present tol is zero ) */
/*        if ( restart ) generate a new starting vector. */
/*     2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}]; */
/*        p_{j} = p_{j}/betaj */
/*     3) r_{j} = OP*v_{j} where OP is defined as in dnaupd */
/*        For shift-invert mode p_{j} = B*v_{j} is already available. */
/*        wnorm = || OP*v_{j} || */
/*     4) Compute the j-th step residual vector. */
/*        w_{j} =  V_{j}^T * B * OP * v_{j} */
/*        r_{j} =  OP*v_{j} - V_{j} * w_{j} */
/*        H(:,j) = w_{j}; */
/*        H(j,j-1) = rnorm */
/*        rnorm = || r_(j) || */
/*        If (rnorm > 0.717*wnorm) accept step and go back to 1) */
/*     5) Re-orthogonalization step: */
/*        s = V_{j}'*B*r_{j} */
/*        r_{j} = r_{j} - V_{j}*s;  rnorm1 = || r_{j} || */
/*        alphaj = alphaj + s_{j}; */
/*     6) Iterative refinement step: */
/*        If (rnorm1 > 0.717*rnorm) then */
/*           rnorm = rnorm1 */
/*           accept step and go back to 1) */
/*        Else */
/*           rnorm = rnorm1 */
/*           If this is the first time in step 6), go to 5) */
/*           Else r_{j} lies in the span of V_{j} numerically. */
/*              Set r_{j} = 0 and rnorm = 0; go to 1) */
/*        EndIf */
/*  End Do */

/* \EndLib */

/* ----------------------------------------------------------------------- */

static int
dnaitr_ (integer * ido, char *bmat, integer * n, integer * k,
			integer * np, integer * nb, doublereal * resid, doublereal * rnorm,
			doublereal * v, integer * ldv, doublereal * h__, integer * ldh, integer * ipntr, doublereal * workd, integer * info, ftnlen bmat_len)
{
	/* Initialized data */

	static logical first = TRUE_;

	/* System generated locals */
	integer h_dim1, h_offset, v_dim1, v_offset, i__1, i__2;
	doublereal d__1, d__2;

	/* Builtin functions */
	double sqrt (doublereal);

	/* Local variables */
	static integer i__, j;
	static real t0, t1, t2, t3, t4, t5;
	static integer jj, ipj, irj, ivj;
	static doublereal ulp, tst1;
	static integer ierr, iter;
	static doublereal unfl, ovfl;
	static integer itry;
	static doublereal temp1;
	static logical orth1, orth2, step3, step4;
	static doublereal betaj;
	static integer infol;
	static doublereal xtemp[2];
	static doublereal wnorm;
	static doublereal rnorm1;
	static logical rstart;
	static integer msglvl;
	static doublereal smlnum;


	/*     %----------------------------------------------------% */
	/*     | Include files for debugging and timing information | */
	/*     %----------------------------------------------------% */


	/* \SCCS Information: @(#) */
	/* FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 */

	/*     %---------------------------------% */
	/*     | See debug.doc for documentation | */
	/*     %---------------------------------% */

	/*     %------------------% */
	/*     | Scalar Arguments | */
	/*     %------------------% */

	/*     %--------------------------------% */
	/*     | See stat.doc for documentation | */
	/*     %--------------------------------% */

	/* \SCCS Information: @(#) */
	/* FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 */



	/*     %-----------------% */
	/*     | Array Arguments | */
	/*     %-----------------% */


	/*     %------------% */
	/*     | Parameters | */
	/*     %------------% */


	/*     %---------------% */
	/*     | Local Scalars | */
	/*     %---------------% */


	/*     %-----------------------% */
	/*     | Local Array Arguments | */
	/*     %-----------------------% */


	/*     %----------------------% */
	/*     | External Subroutines | */
	/*     %----------------------% */


	/*     %--------------------% */
	/*     | External Functions | */
	/*     %--------------------% */


	/*     %---------------------% */
	/*     | Intrinsic Functions | */
	/*     %---------------------% */


	/*     %-----------------% */
	/*     | Data statements | */
	/*     %-----------------% */

	/* Parameter adjustments */
	--workd;
	--resid;
	v_dim1 = *ldv;
	v_offset = 1 + v_dim1;
	v -= v_offset;
	h_dim1 = *ldh;
	h_offset = 1 + h_dim1;
	h__ -= h_offset;
	--ipntr;

	/* Function Body */

	/*     %-----------------------% */
	/*     | Executable Statements | */
	/*     %-----------------------% */

	if (first)
	{

		/*        %-----------------------------------------% */
		/*        | Set machine-dependent constants for the | */
		/*        | the splitting and deflation criterion.  | */
		/*        | If norm(H) <= sqrt(OVFL),               | */
		/*        | overflow should not occur.              | */
		/*        | REFERENCE: LAPACK subroutine dlahqr     | */
		/*        %-----------------------------------------% */

		unfl = dlamch_ ("safe minimum", (ftnlen) 12);
		ovfl = 1. / unfl;
		dlabad_ (&unfl, &ovfl);
		ulp = dlamch_ ("precision", (ftnlen) 9);
		smlnum = unfl * (*n / ulp);
		first = FALSE_;
	}

	if (*ido == 0)
	{

		/*        %-------------------------------% */
		/*        | Initialize timing statistics  | */
		/*        | & message level for debugging | */
		/*        %-------------------------------% */

		second_ (&t0);
		msglvl = debug_1.mnaitr;

		/*        %------------------------------% */
		/*        | Initial call to this routine | */
		/*        %------------------------------% */

		*info = 0;
		step3 = FALSE_;
		step4 = FALSE_;
		rstart = FALSE_;
		orth1 = FALSE_;
		orth2 = FALSE_;
		j = *k + 1;
		ipj = 1;
		irj = ipj + *n;
		ivj = irj + *n;
	}

	/*     %-------------------------------------------------% */
	/*     | When in reverse communication mode one of:      | */
	/*     | STEP3, STEP4, ORTH1, ORTH2, RSTART              | */
	/*     | will be .true. when ....                        | */
	/*     | STEP3: return from computing OP*v_{j}.          | */
	/*     | STEP4: return from computing B-norm of OP*v_{j} | */
	/*     | ORTH1: return from computing B-norm of r_{j+1}  | */
	/*     | ORTH2: return from computing B-norm of          | */
	/*     |        correction to the residual vector.       | */
	/*     | RSTART: return from OP computations needed by   | */
	/*     |         dgetv0.                                 | */
	/*     %-------------------------------------------------% */

	if (step3)
	{
		goto L50;
	}
	if (step4)
	{
		goto L60;
	}
	if (orth1)
	{
		goto L70;
	}
	if (orth2)
	{
		goto L90;
	}
	if (rstart)
	{
		goto L30;
	}

	/*     %-----------------------------% */
	/*     | Else this is the first step | */
	/*     %-----------------------------% */

	/*     %--------------------------------------------------------------% */
	/*     |                                                              | */
	/*     |        A R N O L D I     I T E R A T I O N     L O O P       | */
	/*     |                                                              | */
	/*     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) | */
	/*     %--------------------------------------------------------------% */
 L1000:

	if (msglvl > 1)
	{
		ivout_ (&debug_1.logfil, &c__1, &j, &debug_1.ndigit, "_naitr: generat" "ing Arnoldi vector number", (ftnlen) 40);
		dvout_ (&debug_1.logfil, &c__1, rnorm, &debug_1.ndigit, "_naitr: B-no" "rm of the current residual is", (ftnlen) 41);
	}

	/*        %---------------------------------------------------% */
	/*        | STEP 1: Check if the B norm of j-th residual      | */
	/*        | vector is zero. Equivalent to determing whether   | */
	/*        | an exact j-step Arnoldi factorization is present. | */
	/*        %---------------------------------------------------% */

	betaj = *rnorm;
	if (*rnorm > 0.)
	{
		goto L40;
	}

	/*           %---------------------------------------------------% */
	/*           | Invariant subspace found, generate a new starting | */
	/*           | vector which is orthogonal to the current Arnoldi | */
	/*           | basis and continue the iteration.                 | */
	/*           %---------------------------------------------------% */

	if (msglvl > 0)
	{
		ivout_ (&debug_1.logfil, &c__1, &j, &debug_1.ndigit, "_naitr: ****** " "RESTART AT STEP ******", (ftnlen) 37);
	}

	/*           %---------------------------------------------% */
	/*           | ITRY is the loop variable that controls the | */
	/*           | maximum amount of times that a restart is   | */
	/*           | attempted. NRSTRT is used by stat.h         | */
	/*           %---------------------------------------------% */

	betaj = 0.;
	++timing_1.nrstrt;
	itry = 1;
 L20:
	rstart = TRUE_;
	*ido = 0;
 L30:

	/*           %--------------------------------------% */
	/*           | If in reverse communication mode and | */
	/*           | RSTART = .true. flow returns here.   | */
	/*           %--------------------------------------% */

	dgetv0_ (ido, bmat, &itry, &c_false, n, &j, &v[v_offset], ldv, &resid[1], rnorm, &ipntr[1], &workd[1], &ierr, (ftnlen) 1);
	if (*ido != 99)
	{
		goto L9000;
	}
	if (ierr < 0)
	{
		++itry;
		if (itry <= 3)
		{
			goto L20;
		}

		/*              %------------------------------------------------% */
		/*              | Give up after several restart attempts.        | */
		/*              | Set INFO to the size of the invariant subspace | */
		/*              | which spans OP and exit.                       | */
		/*              %------------------------------------------------% */

		*info = j - 1;
		second_ (&t1);
		timing_1.tnaitr += t1 - t0;
		*ido = 99;
		goto L9000;
	}

 L40:

	/*        %---------------------------------------------------------% */
	/*        | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  | */
	/*        | Note that p_{j} = B*r_{j-1}. In order to avoid overflow | */
	/*        | when reciprocating a small RNORM, test against lower    | */
	/*        | machine bound.                                          | */
	/*        %---------------------------------------------------------% */

	dcopy_ (n, &resid[1], &c__1, &v[j * v_dim1 + 1], &c__1);
	if (*rnorm >= unfl)
	{
		temp1 = 1. / *rnorm;
		dscal_ (n, &temp1, &v[j * v_dim1 + 1], &c__1);
		dscal_ (n, &temp1, &workd[ipj], &c__1);
	}
	else
	{

		/*            %-----------------------------------------% */
		/*            | To scale both v_{j} and p_{j} carefully | */
		/*            | use LAPACK routine SLASCL               | */
		/*            %-----------------------------------------% */

		dlascl_ ("General", &i__, &i__, rnorm, &c_b348, n, &c__1, &v[j * v_dim1 + 1], n, &infol, (ftnlen) 7);
		dlascl_ ("General", &i__, &i__, rnorm, &c_b348, n, &c__1, &workd[ipj], n, &infol, (ftnlen) 7);
	}

	/*        %------------------------------------------------------% */
	/*        | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} | */
	/*        | Note that this is not quite yet r_{j}. See STEP 4    | */
	/*        %------------------------------------------------------% */

	step3 = TRUE_;
	++timing_1.nopx;
	second_ (&t2);
	dcopy_ (n, &v[j * v_dim1 + 1], &c__1, &workd[ivj], &c__1);
	ipntr[1] = ivj;
	ipntr[2] = irj;
	ipntr[3] = ipj;
	*ido = 1;

	/*        %-----------------------------------% */
	/*        | Exit in order to compute OP*v_{j} | */
	/*        %-----------------------------------% */

	goto L9000;
 L50:

	/*        %----------------------------------% */
	/*        | Back from reverse communication; | */
	/*        | WORKD(IRJ:IRJ+N-1) := OP*v_{j}   | */
	/*        | if step3 = .true.                | */
	/*        %----------------------------------% */

	second_ (&t3);
	timing_1.tmvopx += t3 - t2;
	step3 = FALSE_;

	/*        %------------------------------------------% */
	/*        | Put another copy of OP*v_{j} into RESID. | */
	/*        %------------------------------------------% */

	dcopy_ (n, &workd[irj], &c__1, &resid[1], &c__1);

	/*        %---------------------------------------% */
	/*        | STEP 4:  Finish extending the Arnoldi | */
	/*        |          factorization to length j.   | */
	/*        %---------------------------------------% */

	second_ (&t2);
	if (*(unsigned char *) bmat == 'G')
	{
		++timing_1.nbx;
		step4 = TRUE_;
		ipntr[1] = irj;
		ipntr[2] = ipj;
		*ido = 2;

		/*           %-------------------------------------% */
		/*           | Exit in order to compute B*OP*v_{j} | */
		/*           %-------------------------------------% */

		goto L9000;
	}
	else if (*(unsigned char *) bmat == 'I')
	{
		dcopy_ (n, &resid[1], &c__1, &workd[ipj], &c__1);
	}
 L60:

	/*        %----------------------------------% */
	/*        | Back from reverse communication; | */
	/*        | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j} | */
	/*        | if step4 = .true.                | */
	/*        %----------------------------------% */

	if (*(unsigned char *) bmat == 'G')
	{
		second_ (&t3);
		timing_1.tmvbx += t3 - t2;
	}

	step4 = FALSE_;

	/*        %-------------------------------------% */
	/*        | The following is needed for STEP 5. | */
	/*        | Compute the B-norm of OP*v_{j}.     | */
	/*        %-------------------------------------% */

	if (*(unsigned char *) bmat == 'G')
	{
		wnorm = ddot_ (n, &resid[1], &c__1, &workd[ipj], &c__1);
		wnorm = sqrt ((abs (wnorm)));
	}
	else if (*(unsigned char *) bmat == 'I')
	{
		wnorm = dnrm2_ (n, &resid[1], &c__1);
	}

	/*        %-----------------------------------------% */
	/*        | Compute the j-th residual corresponding | */
	/*        | to the j step factorization.            | */
	/*        | Use Classical Gram Schmidt and compute: | */
	/*        | w_{j} <-  V_{j}^T * B * OP * v_{j}      | */
	/*        | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      | */
	/*        %-----------------------------------------% */


	/*        %------------------------------------------% */
	/*        | Compute the j Fourier coefficients w_{j} | */
	/*        | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  | */
	/*        %------------------------------------------% */

	dgemv_ ("T", n, &j, &c_b348, &v[v_offset], ldv, &workd[ipj], &c__1, &c_b507, &h__[j * h_dim1 + 1], &c__1, (ftnlen) 1);

	/*        %--------------------------------------% */
	/*        | Orthogonalize r_{j} against V_{j}.   | */
	/*        | RESID contains OP*v_{j}. See STEP 3. | */
	/*        %--------------------------------------% */

	dgemv_ ("N", n, &j, &c_b347, &v[v_offset], ldv, &h__[j * h_dim1 + 1], &c__1, &c_b348, &resid[1], &c__1, (ftnlen) 1);

	if (j > 1)
	{
		h__[j + (j - 1) * h_dim1] = betaj;
	}

	second_ (&t4);

	orth1 = TRUE_;

	second_ (&t2);
	if (*(unsigned char *) bmat == 'G')
	{
		++timing_1.nbx;
		dcopy_ (n, &resid[1], &c__1, &workd[irj], &c__1);
		ipntr[1] = irj;
		ipntr[2] = ipj;
		*ido = 2;

		/*           %----------------------------------% */
		/*           | Exit in order to compute B*r_{j} | */
		/*           %----------------------------------% */

		goto L9000;
	}
	else if (*(unsigned char *) bmat == 'I')
	{
		dcopy_ (n, &resid[1], &c__1, &workd[ipj], &c__1);
	}
 L70:

	/*        %---------------------------------------------------% */
	/*        | Back from reverse communication if ORTH1 = .true. | */
	/*        | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    | */
	/*        %---------------------------------------------------% */

	if (*(unsigned char *) bmat == 'G')
	{
		second_ (&t3);
		timing_1.tmvbx += t3 - t2;
	}

	orth1 = FALSE_;

	/*        %------------------------------% */
	/*        | Compute the B-norm of r_{j}. | */
	/*        %------------------------------% */

	if (*(unsigned char *) bmat == 'G')
	{
		*rnorm = ddot_ (n, &resid[1], &c__1, &workd[ipj], &c__1);
		*rnorm = sqrt ((abs (*rnorm)));
	}
	else if (*(unsigned char *) bmat == 'I')
	{
		*rnorm = dnrm2_ (n, &resid[1], &c__1);
	}

	/*        %-----------------------------------------------------------% */
	/*        | STEP 5: Re-orthogonalization / Iterative refinement phase | */
	/*        | Maximum NITER_ITREF tries.                                | */
	/*        |                                                           | */
	/*        |          s      = V_{j}^T * B * r_{j}                     | */
	/*        |          r_{j}  = r_{j} - V_{j}*s                         | */
	/*        |          alphaj = alphaj + s_{j}                          | */
	/*        |                                                           | */
	/*        | The stopping criteria used for iterative refinement is    | */
	/*        | discussed in Parlett's book SEP, page 107 and in Gragg &  | */
	/*        | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         | */
	/*        | Determine if we need to correct the residual. The goal is | */
	/*        | to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  | */
	/*        | The following test determines whether the sine of the     | */
	/*        | angle between  OP*x and the computed residual is less     | */
	/*        | than or equal to 0.717.                                   | */
	/*        %-----------------------------------------------------------% */

	if (*rnorm > wnorm * .717f)
	{
		goto L100;
	}
	iter = 0;
	++timing_1.nrorth;

	/*        %---------------------------------------------------% */
	/*        | Enter the Iterative refinement phase. If further  | */
	/*        | refinement is necessary, loop back here. The loop | */
	/*        | variable is ITER. Perform a step of Classical     | */
	/*        | Gram-Schmidt using all the Arnoldi vectors V_{j}  | */
	/*        %---------------------------------------------------% */

 L80:

	if (msglvl > 2)
	{
		xtemp[0] = wnorm;
		xtemp[1] = *rnorm;
		dvout_ (&debug_1.logfil, &c__2, xtemp, &debug_1.ndigit, "_naitr: re-o" "rthonalization; wnorm and rnorm are", (ftnlen) 47);
		dvout_ (&debug_1.logfil, &j, &h__[j * h_dim1 + 1], &debug_1.ndigit, "_naitr: j-th column of H", (ftnlen) 24);
	}

	/*        %----------------------------------------------------% */
	/*        | Compute V_{j}^T * B * r_{j}.                       | */
	/*        | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). | */
	/*        %----------------------------------------------------% */

	dgemv_ ("T", n, &j, &c_b348, &v[v_offset], ldv, &workd[ipj], &c__1, &c_b507, &workd[irj], &c__1, (ftnlen) 1);

	/*        %---------------------------------------------% */
	/*        | Compute the correction to the residual:     | */
	/*        | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1). | */
	/*        | The correction to H is v(:,1:J)*H(1:J,1:J)  | */
	/*        | + v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j.         | */
	/*        %---------------------------------------------% */

	dgemv_ ("N", n, &j, &c_b347, &v[v_offset], ldv, &workd[irj], &c__1, &c_b348, &resid[1], &c__1, (ftnlen) 1);
	daxpy_ (&j, &c_b348, &workd[irj], &c__1, &h__[j * h_dim1 + 1], &c__1);

	orth2 = TRUE_;
	second_ (&t2);
	if (*(unsigned char *) bmat == 'G')
	{
		++timing_1.nbx;
		dcopy_ (n, &resid[1], &c__1, &workd[irj], &c__1);
		ipntr[1] = irj;
		ipntr[2] = ipj;
		*ido = 2;

		/*           %-----------------------------------% */
		/*           | Exit in order to compute B*r_{j}. | */
		/*           | r_{j} is the corrected residual.  | */
		/*           %-----------------------------------% */

		goto L9000;
	}
	else if (*(unsigned char *) bmat == 'I')
	{
		dcopy_ (n, &resid[1], &c__1, &workd[ipj], &c__1);
	}
 L90:

	/*        %---------------------------------------------------% */
	/*        | Back from reverse communication if ORTH2 = .true. | */
	/*        %---------------------------------------------------% */

	if (*(unsigned char *) bmat == 'G')
	{
		second_ (&t3);
		timing_1.tmvbx += t3 - t2;
	}

	/*        %-----------------------------------------------------% */
	/*        | Compute the B-norm of the corrected residual r_{j}. | */
	/*        %-----------------------------------------------------% */

	if (*(unsigned char *) bmat == 'G')
	{
		rnorm1 = ddot_ (n, &resid[1], &c__1, &workd[ipj], &c__1);
		rnorm1 = sqrt ((abs (rnorm1)));
	}
	else if (*(unsigned char *) bmat == 'I')
	{
		rnorm1 = dnrm2_ (n, &resid[1], &c__1);
	}

	if (msglvl > 0 && iter > 0)
	{
		ivout_ (&debug_1.logfil, &c__1, &j, &debug_1.ndigit, "_naitr: Iterati" "ve refinement for Arnoldi residual", (ftnlen) 49);
		if (msglvl > 2)
		{
			xtemp[0] = *rnorm;
			xtemp[1] = rnorm1;
			dvout_ (&debug_1.logfil, &c__2, xtemp, &debug_1.ndigit, "_naitr: " "iterative refinement ; rnorm and rnorm1 are", (ftnlen) 51);
		}
	}

	/*        %-----------------------------------------% */
	/*        | Determine if we need to perform another | */
	/*        | step of re-orthogonalization.           | */
	/*        %-----------------------------------------% */

	if (rnorm1 > *rnorm * .717f)
	{

		/*           %---------------------------------------% */
		/*           | No need for further refinement.       | */
		/*           | The cosine of the angle between the   | */
		/*           | corrected residual vector and the old | */
		/*           | residual vector is greater than 0.717 | */
		/*           | In other words the corrected residual | */
		/*           | and the old residual vector share an  | */
		/*           | angle of less than arcCOS(0.717)      | */
		/*           %---------------------------------------% */

		*rnorm = rnorm1;

	}
	else
	{

		/*           %-------------------------------------------% */
		/*           | Another step of iterative refinement step | */
		/*           | is required. NITREF is used by stat.h     | */
		/*           %-------------------------------------------% */

		++timing_1.nitref;
		*rnorm = rnorm1;
		++iter;
		if (iter <= 1)
		{
			goto L80;
		}

		/*           %-------------------------------------------------% */
		/*           | Otherwise RESID is numerically in the span of V | */
		/*           %-------------------------------------------------% */

		i__1 = *n;
		for (jj = 1; jj <= i__1; ++jj)
		{
			resid[jj] = 0.;
			/* L95: */
		}
		*rnorm = 0.;
	}

	/*        %----------------------------------------------% */
	/*        | Branch here directly if iterative refinement | */
	/*        | wasn't necessary or after at most NITER_REF  | */
	/*        | steps of iterative refinement.               | */
	/*        %----------------------------------------------% */

 L100:

	rstart = FALSE_;
	orth2 = FALSE_;

	second_ (&t5);
	timing_1.titref += t5 - t4;

	/*        %------------------------------------% */
	/*        | STEP 6: Update  j = j+1;  Continue | */
	/*        %------------------------------------% */

	++j;
	if (j > *k + *np)
	{
		second_ (&t1);
		timing_1.tnaitr += t1 - t0;
		*ido = 99;
		i__1 = *k + *np - 1;
		for (i__ = max (1, *k); i__ <= i__1; ++i__)
		{

			/*              %--------------------------------------------% */
			/*              | Check for splitting and deflation.         | */
			/*              | Use a standard test as in the QR algorithm | */
			/*              | REFERENCE: LAPACK subroutine dlahqr        | */
			/*              %--------------------------------------------% */

			tst1 = (d__1 = h__[i__ + i__ * h_dim1], abs (d__1)) + (d__2 = h__[i__ + 1 + (i__ + 1) * h_dim1], abs (d__2));
			if (tst1 == 0.)
			{
				i__2 = *k + *np;
				tst1 = dlanhs_ ("1", &i__2, &h__[h_offset], ldh, &workd[*n + 1], (ftnlen) 1);
			}
			/* Computing MAX */
			d__2 = ulp * tst1;
			if ((d__1 = h__[i__ + 1 + i__ * h_dim1], abs (d__1)) <= max (d__2, smlnum))
			{
				h__[i__ + 1 + i__ * h_dim1] = 0.;
			}
			/* L110: */
		}

		if (msglvl > 2)
		{
			i__1 = *k + *np;
			i__2 = *k + *np;
			dmout_ (&debug_1.logfil, &i__1, &i__2, &h__[h_offset], ldh, &debug_1.ndigit, "_naitr: Final upper Hessenberg matrix H" " of order K+NP", (ftnlen) 53);
		}

		goto L9000;
	}

	/*        %--------------------------------------------------------% */
	/*        | Loop back to extend the factorization by another step. | */
	/*        %--------------------------------------------------------% */

	goto L1000;

	/*     %---------------------------------------------------------------% */
	/*     |                                                               | */
	/*     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  | */
	/*     |                                                               | */
	/*     %---------------------------------------------------------------% */

 L9000:
	return 0;

	/*     %---------------% */
	/*     | End of dnaitr | */
	/*     %---------------% */

}	/* dnaitr_ */

/* ----------------------------------------------------------------------- */
/* \BeginDoc */

/* \Name: dnapps */

/* \Description: */
/*  Given the Arnoldi factorization */

/*     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T, */

/*  apply NP implicit shifts resulting in */

/*     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q */

/*  where Q is an orthogonal matrix which is the product of rotations */
/*  and reflections resulting from the NP bulge chage sweeps. */
/*  The updated Arnoldi factorization becomes: */

/*     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T. */

/* \Usage: */
/*  call dnapps */
/*     ( N, KEV, NP, SHIFTR, SHIFTI, V, LDV, H, LDH, RESID, Q, LDQ, */
/*       WORKL, WORKD ) */

/* \Arguments */
/*  N       Integer.  (INPUT) */
/*          Problem size, i.e. size of matrix A. */

/*  KEV     Integer.  (INPUT/OUTPUT) */
/*          KEV+NP is the size of the input matrix H. */
/*          KEV is the size of the updated matrix HNEW.  KEV is only */
/*          updated on ouput when fewer than NP shifts are applied in */
/*          order to keep the conjugate pair together. */

/*  NP      Integer.  (INPUT) */
/*          Number of implicit shifts to be applied. */

/*  SHIFTR, Double precision array of length NP.  (INPUT) */
/*  SHIFTI  Real and imaginary part of the shifts to be applied. */
/*          Upon, entry to dnapps, the shifts must be sorted so that the */
/*          conjugate pairs are in consecutive locations. */

/*  V       Double precision N by (KEV+NP) array.  (INPUT/OUTPUT) */
/*          On INPUT, V contains the current KEV+NP Arnoldi vectors. */
/*          On OUTPUT, V contains the updated KEV Arnoldi vectors */
/*          in the first KEV columns of V. */

/*  LDV     Integer.  (INPUT) */
/*          Leading dimension of V exactly as declared in the calling */
/*          program. */

/*  H       Double precision (KEV+NP) by (KEV+NP) array.  (INPUT/OUTPUT) */
/*          On INPUT, H contains the current KEV+NP by KEV+NP upper */
/*          Hessenber matrix of the Arnoldi factorization. */
/*          On OUTPUT, H contains the updated KEV by KEV upper Hessenberg */
/*          matrix in the KEV leading submatrix. */

/*  LDH     Integer.  (INPUT) */
/*          Leading dimension of H exactly as declared in the calling */
/*          program. */

/*  RESID   Double precision array of length N.  (INPUT/OUTPUT) */
/*          On INPUT, RESID contains the the residual vector r_{k+p}. */
/*          On OUTPUT, RESID is the update residual vector rnew_{k} */
/*          in the first KEV locations. */

/*  Q       Double precision KEV+NP by KEV+NP work array.  (WORKSPACE) */
/*          Work array used to accumulate the rotations and reflections */
/*          during the bulge chase sweep. */

/*  LDQ     Integer.  (INPUT) */
/*          Leading dimension of Q exactly as declared in the calling */
/*          program. */

/*  WORKL   Double precision work array of length (KEV+NP).  (WORKSPACE) */
/*          Private (replicated) array on each PE or array allocated on */
/*          the front end. */

/*  WORKD   Double precision work array of length 2*N.  (WORKSPACE) */
/*          Distributed array used in the application of the accumulated */
/*          orthogonal matrix Q. */

/* \EndDoc */

/* ----------------------------------------------------------------------- */

/* \BeginLib */

/* \Local variables: */
/*     xxxxxx  real */

/* \References: */
/*  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in */
/*     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992), */
/*     pp 357-385. */

/* \Routines called: */
/*     ivout   ARPACK utility routine that prints integers. */
/*     second  ARPACK utility routine for timing. */
/*     dmout   ARPACK utility routine that prints matrices. */
/*     dvout   ARPACK utility routine that prints vectors. */
/*     dlabad  LAPACK routine that computes machine constants. */
/*     dlacpy  LAPACK matrix copy routine. */
/*     dlamch  LAPACK routine that determines machine constants. */
/*     dlanhs  LAPACK routine that computes various norms of a matrix. */
/*     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully. */
/*     dlarf   LAPACK routine that applies Householder reflection to */
/*             a matrix. */
/*     dlarfg  LAPACK Householder reflection construction routine. */
/*     dlartg  LAPACK Givens rotation construction routine. */
/*     dlaset  LAPACK matrix initialization routine. */
/*     dgemv   Level 2 BLAS routine for matrix vector multiplication. */
/*     daxpy   Level 1 BLAS that computes a vector triad. */
/*     dcopy   Level 1 BLAS that copies one vector to another . */
/*     dscal   Level 1 BLAS that scales a vector. */

/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \Revision history: */
/*     xx/xx/92: Version ' 2.4' */

/* \SCCS Information: @(#) */
/* FILE: napps.F   SID: 2.4   DATE OF SID: 3/28/97   RELEASE: 2 */

/* \Remarks */
/*  1. In this version, each shift is applied to all the sublocks of */
/*     the Hessenberg matrix H and not just to the submatrix that it */
/*     comes from. Deflation as in LAPACK routine dlahqr (QR algorithm */
/*     for upper Hessenberg matrices ) is used. */
/*     The subdiagonals of H are enforced to be non-negative. */

/* \EndLib */

/* ----------------------------------------------------------------------- */

static int
dnapps_ (integer * n, integer * kev, integer * np,
			doublereal * shiftr, doublereal * shifti, doublereal * v, integer * ldv,
			doublereal * h__, integer * ldh, doublereal * resid, doublereal * q, integer * ldq, doublereal * workl, doublereal * workd)
{
	/* Initialized data */

	static logical first = TRUE_;

	/* System generated locals */
	integer h_dim1, h_offset, v_dim1, v_offset, q_dim1, q_offset, i__1, i__2, i__3, i__4;
	doublereal d__1, d__2;

	/* Local variables */
	static doublereal c__, f, g;
	static integer i__, j;
	static doublereal r__, s, t, u[3];
	static real t0, t1;
	static doublereal h11, h12, h21, h22, h32;
	static integer jj, ir, nr;
	static doublereal tau, ulp, tst1;
	static integer iend;
	static doublereal unfl, ovfl;
	static logical cconj;
	static doublereal sigmai;
	static integer istart, kplusp, msglvl;
	static doublereal sigmar, smlnum;


	/*     %----------------------------------------------------% */
	/*     | Include files for debugging and timing information | */
	/*     %----------------------------------------------------% */


	/* \SCCS Information: @(#) */
	/* FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 */

	/*     %---------------------------------% */
	/*     | See debug.doc for documentation | */
	/*     %---------------------------------% */

	/*     %------------------% */
	/*     | Scalar Arguments | */
	/*     %------------------% */

	/*     %--------------------------------% */
	/*     | See stat.doc for documentation | */
	/*     %--------------------------------% */

	/* \SCCS Information: @(#) */
	/* FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 */



	/*     %-----------------% */
	/*     | Array Arguments | */
	/*     %-----------------% */


	/*     %------------% */
	/*     | Parameters | */
	/*     %------------% */


	/*     %------------------------% */
	/*     | Local Scalars & Arrays | */
	/*     %------------------------% */


	/*     %----------------------% */
	/*     | External Subroutines | */
	/*     %----------------------% */


	/*     %--------------------% */
	/*     | External Functions | */
	/*     %--------------------% */


	/*     %----------------------% */
	/*     | Intrinsics Functions | */
	/*     %----------------------% */


	/*     %----------------% */
	/*     | Data statments | */
	/*     %----------------% */

	/* Parameter adjustments */
	--workd;
	--resid;
	--workl;
	--shifti;
	--shiftr;
	v_dim1 = *ldv;
	v_offset = 1 + v_dim1;
	v -= v_offset;
	h_dim1 = *ldh;
	h_offset = 1 + h_dim1;
	h__ -= h_offset;
	q_dim1 = *ldq;
	q_offset = 1 + q_dim1;
	q -= q_offset;

	/* Function Body */

	/*     %-----------------------% */
	/*     | Executable Statements | */
	/*     %-----------------------% */

	if (first)
	{

		/*        %-----------------------------------------------% */
		/*        | Set machine-dependent constants for the       | */
		/*        | stopping criterion. If norm(H) <= sqrt(OVFL), | */
		/*        | overflow should not occur.                    | */
		/*        | REFERENCE: LAPACK subroutine dlahqr           | */
		/*        %-----------------------------------------------% */

		unfl = dlamch_ ("safe minimum", (ftnlen) 12);
		ovfl = 1. / unfl;
		dlabad_ (&unfl, &ovfl);
		ulp = dlamch_ ("precision", (ftnlen) 9);
		smlnum = unfl * (*n / ulp);
		first = FALSE_;
	}

	/*     %-------------------------------% */
	/*     | Initialize timing statistics  | */
	/*     | & message level for debugging | */
	/*     %-------------------------------% */

	second_ (&t0);
	msglvl = debug_1.mnapps;
	kplusp = *kev + *np;

	/*     %--------------------------------------------% */
	/*     | Initialize Q to the identity to accumulate | */
	/*     | the rotations and reflections              | */
	/*     %--------------------------------------------% */

	dlaset_ ("All", &kplusp, &kplusp, &c_b507, &c_b348, &q[q_offset], ldq, (ftnlen) 3);

	/*     %----------------------------------------------% */
	/*     | Quick return if there are no shifts to apply | */
	/*     %----------------------------------------------% */

	if (*np == 0)
	{
		goto L9000;
	}

	/*     %----------------------------------------------% */
	/*     | Chase the bulge with the application of each | */
	/*     | implicit shift. Each shift is applied to the | */
	/*     | whole matrix including each block.           | */
	/*     %----------------------------------------------% */

	cconj = FALSE_;
	i__1 = *np;
	for (jj = 1; jj <= i__1; ++jj)
	{
		sigmar = shiftr[jj];
		sigmai = shifti[jj];

		if (msglvl > 2)
		{
			ivout_ (&debug_1.logfil, &c__1, &jj, &debug_1.ndigit, "_napps: sh" "ift number.", (ftnlen) 21);
			dvout_ (&debug_1.logfil, &c__1, &sigmar, &debug_1.ndigit, "_napps" ": The real part of the shift ", (ftnlen) 35);
			dvout_ (&debug_1.logfil, &c__1, &sigmai, &debug_1.ndigit, "_napps" ": The imaginary part of the shift ", (ftnlen) 40);
		}

		/*        %-------------------------------------------------% */
		/*        | The following set of conditionals is necessary  | */
		/*        | in order that complex conjugate pairs of shifts | */
		/*        | are applied together or not at all.             | */
		/*        %-------------------------------------------------% */

		if (cconj)
		{

			/*           %-----------------------------------------% */
			/*           | cconj = .true. means the previous shift | */
			/*           | had non-zero imaginary part.            | */
			/*           %-----------------------------------------% */

			cconj = FALSE_;
			goto L110;
		}
		else if (jj < *np && abs (sigmai) > 0.)
		{

			/*           %------------------------------------% */
			/*           | Start of a complex conjugate pair. | */
			/*           %------------------------------------% */

			cconj = TRUE_;
		}
		else if (jj == *np && abs (sigmai) > 0.)
		{

			/*           %----------------------------------------------% */
			/*           | The last shift has a nonzero imaginary part. | */
			/*           | Don't apply it; thus the order of the        | */
			/*           | compressed H is order KEV+1 since only np-1  | */
			/*           | were applied.                                | */
			/*           %----------------------------------------------% */

			++(*kev);
			goto L110;
		}
		istart = 1;
	 L20:

		/*        %--------------------------------------------------% */
		/*        | if sigmai = 0 then                               | */
		/*        |    Apply the jj-th shift ...                     | */
		/*        | else                                             | */
		/*        |    Apply the jj-th and (jj+1)-th together ...    | */
		/*        |    (Note that jj < np at this point in the code) | */
		/*        | end                                              | */
		/*        | to the current block of H. The next do loop      | */
		/*        | determines the current block ;                   | */
		/*        %--------------------------------------------------% */

		i__2 = kplusp - 1;
		for (i__ = istart; i__ <= i__2; ++i__)
		{

			/*           %----------------------------------------% */
			/*           | Check for splitting and deflation. Use | */
			/*           | a standard test as in the QR algorithm | */
			/*           | REFERENCE: LAPACK subroutine dlahqr    | */
			/*           %----------------------------------------% */

			tst1 = (d__1 = h__[i__ + i__ * h_dim1], abs (d__1)) + (d__2 = h__[i__ + 1 + (i__ + 1) * h_dim1], abs (d__2));
			if (tst1 == 0.)
			{
				i__3 = kplusp - jj + 1;
				tst1 = dlanhs_ ("1", &i__3, &h__[h_offset], ldh, &workl[1], (ftnlen) 1);
			}
			/* Computing MAX */
			d__2 = ulp * tst1;
			if ((d__1 = h__[i__ + 1 + i__ * h_dim1], abs (d__1)) <= max (d__2, smlnum))
			{
				if (msglvl > 0)
				{
					ivout_ (&debug_1.logfil, &c__1, &i__, &debug_1.ndigit, "_napps: matrix splitting at row/column no.", (ftnlen) 42);
					ivout_ (&debug_1.logfil, &c__1, &jj, &debug_1.ndigit, "_napps: matrix splitting with shift number.", (ftnlen) 43);
					dvout_ (&debug_1.logfil, &c__1, &h__[i__ + 1 + i__ * h_dim1], &debug_1.ndigit, "_napps: off diagonal " "element.", (ftnlen) 29);
				}
				iend = i__;
				h__[i__ + 1 + i__ * h_dim1] = 0.;
				goto L40;
			}
			/* L30: */
		}
		iend = kplusp;
	 L40:

		if (msglvl > 2)
		{
			ivout_ (&debug_1.logfil, &c__1, &istart, &debug_1.ndigit, "_napps" ": Start of current block ", (ftnlen) 31);
			ivout_ (&debug_1.logfil, &c__1, &iend, &debug_1.ndigit, "_napps: " "End of current block ", (ftnlen) 29);
		}

		/*        %------------------------------------------------% */
		/*        | No reason to apply a shift to block of order 1 | */
		/*        %------------------------------------------------% */

		if (istart == iend)
		{
			goto L100;
		}

		/*        %------------------------------------------------------% */
		/*        | If istart + 1 = iend then no reason to apply a       | */
		/*        | complex conjugate pair of shifts on a 2 by 2 matrix. | */
		/*        %------------------------------------------------------% */

		if (istart + 1 == iend && abs (sigmai) > 0.)
		{
			goto L100;
		}

		h11 = h__[istart + istart * h_dim1];
		h21 = h__[istart + 1 + istart * h_dim1];
		if (abs (sigmai) <= 0.)
		{

			/*           %---------------------------------------------% */
			/*           | Real-valued shift ==> apply single shift QR | */
			/*           %---------------------------------------------% */

			f = h11 - sigmar;
			g = h21;

			i__2 = iend - 1;
			for (i__ = istart; i__ <= i__2; ++i__)
			{

				/*              %-----------------------------------------------------% */
				/*              | Contruct the plane rotation G to zero out the bulge | */
				/*              %-----------------------------------------------------% */

				dlartg_ (&f, &g, &c__, &s, &r__);
				if (i__ > istart)
				{

					/*                 %-------------------------------------------% */
					/*                 | The following ensures that h(1:iend-1,1), | */
					/*                 | the first iend-2 off diagonal of elements | */
					/*                 | H, remain non negative.                   | */
					/*                 %-------------------------------------------% */

					if (r__ < 0.)
					{
						r__ = -r__;
						c__ = -c__;
						s = -s;
					}
					h__[i__ + (i__ - 1) * h_dim1] = r__;
					h__[i__ + 1 + (i__ - 1) * h_dim1] = 0.;
				}

				/*              %---------------------------------------------% */
				/*              | Apply rotation to the left of H;  H <- G'*H | */
				/*              %---------------------------------------------% */

				i__3 = kplusp;
				for (j = i__; j <= i__3; ++j)
				{
					t = c__ * h__[i__ + j * h_dim1] + s * h__[i__ + 1 + j * h_dim1];
					h__[i__ + 1 + j * h_dim1] = -s * h__[i__ + j * h_dim1] + c__ * h__[i__ + 1 + j * h_dim1];
					h__[i__ + j * h_dim1] = t;
					/* L50: */
				}

				/*              %---------------------------------------------% */
				/*              | Apply rotation to the right of H;  H <- H*G | */
				/*              %---------------------------------------------% */

				/* Computing MIN */
				i__4 = i__ + 2;
				i__3 = min (i__4, iend);
				for (j = 1; j <= i__3; ++j)
				{
					t = c__ * h__[j + i__ * h_dim1] + s * h__[j + (i__ + 1) * h_dim1];
					h__[j + (i__ + 1) * h_dim1] = -s * h__[j + i__ * h_dim1] + c__ * h__[j + (i__ + 1) * h_dim1];
					h__[j + i__ * h_dim1] = t;
					/* L60: */
				}

				/*              %----------------------------------------------------% */
				/*              | Accumulate the rotation in the matrix Q;  Q <- Q*G | */
				/*              %----------------------------------------------------% */

				/* Computing MIN */
				i__4 = i__ + jj;
				i__3 = min (i__4, kplusp);
				for (j = 1; j <= i__3; ++j)
				{
					t = c__ * q[j + i__ * q_dim1] + s * q[j + (i__ + 1) * q_dim1];
					q[j + (i__ + 1) * q_dim1] = -s * q[j + i__ * q_dim1] + c__ * q[j + (i__ + 1) * q_dim1];
					q[j + i__ * q_dim1] = t;
					/* L70: */
				}

				/*              %---------------------------% */
				/*              | Prepare for next rotation | */
				/*              %---------------------------% */

				if (i__ < iend - 1)
				{
					f = h__[i__ + 1 + i__ * h_dim1];
					g = h__[i__ + 2 + i__ * h_dim1];
				}
				/* L80: */
			}

			/*           %-----------------------------------% */
			/*           | Finished applying the real shift. | */
			/*           %-----------------------------------% */

		}
		else
		{

			/*           %----------------------------------------------------% */
			/*           | Complex conjugate shifts ==> apply double shift QR | */
			/*           %----------------------------------------------------% */

			h12 = h__[istart + (istart + 1) * h_dim1];
			h22 = h__[istart + 1 + (istart + 1) * h_dim1];
			h32 = h__[istart + 2 + (istart + 1) * h_dim1];

			/*           %---------------------------------------------------------% */
			/*           | Compute 1st column of (H - shift*I)*(H - conj(shift)*I) | */
			/*           %---------------------------------------------------------% */

			s = sigmar * 2.f;
			t = dlapy2_ (&sigmar, &sigmai);
			u[0] = (h11 * (h11 - s) + t * t) / h21 + h12;
			u[1] = h11 + h22 - s;
			u[2] = h32;

			i__2 = iend - 1;
			for (i__ = istart; i__ <= i__2; ++i__)
			{

				/* Computing MIN */
				i__3 = 3, i__4 = iend - i__ + 1;
				nr = min (i__3, i__4);

				/*              %-----------------------------------------------------% */
				/*              | Construct Householder reflector G to zero out u(1). | */
				/*              | G is of the form I - tau*( 1 u )' * ( 1 u' ).       | */
				/*              %-----------------------------------------------------% */

				dlarfg_ (&nr, u, &u[1], &c__1, &tau);

				if (i__ > istart)
				{
					h__[i__ + (i__ - 1) * h_dim1] = u[0];
					h__[i__ + 1 + (i__ - 1) * h_dim1] = 0.;
					if (i__ < iend - 1)
					{
						h__[i__ + 2 + (i__ - 1) * h_dim1] = 0.;
					}
				}
				u[0] = 1.;

				/*              %--------------------------------------% */
				/*              | Apply the reflector to the left of H | */
				/*              %--------------------------------------% */

				i__3 = kplusp - i__ + 1;
				dlarf_ ("Left", &nr, &i__3, u, &c__1, &tau, &h__[i__ + i__ * h_dim1], ldh, &workl[1], (ftnlen) 4);

				/*              %---------------------------------------% */
				/*              | Apply the reflector to the right of H | */
				/*              %---------------------------------------% */

				/* Computing MIN */
				i__3 = i__ + 3;
				ir = min (i__3, iend);
				dlarf_ ("Right", &ir, &nr, u, &c__1, &tau, &h__[i__ * h_dim1 + 1], ldh, &workl[1], (ftnlen) 5);

				/*              %-----------------------------------------------------% */
				/*              | Accumulate the reflector in the matrix Q;  Q <- Q*G | */
				/*              %-----------------------------------------------------% */

				dlarf_ ("Right", &kplusp, &nr, u, &c__1, &tau, &q[i__ * q_dim1 + 1], ldq, &workl[1], (ftnlen) 5);

				/*              %----------------------------% */
				/*              | Prepare for next reflector | */
				/*              %----------------------------% */

				if (i__ < iend - 1)
				{
					u[0] = h__[i__ + 1 + i__ * h_dim1];
					u[1] = h__[i__ + 2 + i__ * h_dim1];
					if (i__ < iend - 2)
					{
						u[2] = h__[i__ + 3 + i__ * h_dim1];
					}
				}

				/* L90: */
			}

			/*           %--------------------------------------------% */
			/*           | Finished applying a complex pair of shifts | */
			/*           | to the current block                       | */
			/*           %--------------------------------------------% */

		}

	 L100:

		/*        %---------------------------------------------------------% */
		/*        | Apply the same shift to the next block if there is any. | */
		/*        %---------------------------------------------------------% */

		istart = iend + 1;
		if (iend < kplusp)
		{
			goto L20;
		}

		/*        %---------------------------------------------% */
		/*        | Loop back to the top to get the next shift. | */
		/*        %---------------------------------------------% */

	 L110:
		;
	}

	/*     %--------------------------------------------------% */
	/*     | Perform a similarity transformation that makes   | */
	/*     | sure that H will have non negative sub diagonals | */
	/*     %--------------------------------------------------% */

	i__1 = *kev;
	for (j = 1; j <= i__1; ++j)
	{
		if (h__[j + 1 + j * h_dim1] < 0.)
		{
			i__2 = kplusp - j + 1;
			dscal_ (&i__2, &c_b347, &h__[j + 1 + j * h_dim1], ldh);
			/* Computing MIN */
			i__3 = j + 2;
			i__2 = min (i__3, kplusp);
			dscal_ (&i__2, &c_b347, &h__[(j + 1) * h_dim1 + 1], &c__1);
			/* Computing MIN */
			i__3 = j + *np + 1;
			i__2 = min (i__3, kplusp);
			dscal_ (&i__2, &c_b347, &q[(j + 1) * q_dim1 + 1], &c__1);
		}
		/* L120: */
	}

	i__1 = *kev;
	for (i__ = 1; i__ <= i__1; ++i__)
	{

		/*        %--------------------------------------------% */
		/*        | Final check for splitting and deflation.   | */
		/*        | Use a standard test as in the QR algorithm | */
		/*        | REFERENCE: LAPACK subroutine dlahqr        | */
		/*        %--------------------------------------------% */

		tst1 = (d__1 = h__[i__ + i__ * h_dim1], abs (d__1)) + (d__2 = h__[i__ + 1 + (i__ + 1) * h_dim1], abs (d__2));
		if (tst1 == 0.)
		{
			tst1 = dlanhs_ ("1", kev, &h__[h_offset], ldh, &workl[1], (ftnlen) 1);
		}
		/* Computing MAX */
		d__1 = ulp * tst1;
		if (h__[i__ + 1 + i__ * h_dim1] <= max (d__1, smlnum))
		{
			h__[i__ + 1 + i__ * h_dim1] = 0.;
		}
		/* L130: */
	}

	/*     %-------------------------------------------------% */
	/*     | Compute the (kev+1)-st column of (V*Q) and      | */
	/*     | temporarily store the result in WORKD(N+1:2*N). | */
	/*     | This is needed in the residual update since we  | */
	/*     | cannot GUARANTEE that the corresponding entry   | */
	/*     | of H would be zero as in exact arithmetic.      | */
	/*     %-------------------------------------------------% */

	if (h__[*kev + 1 + *kev * h_dim1] > 0.)
	{
		dgemv_ ("N", n, &kplusp, &c_b348, &v[v_offset], ldv, &q[(*kev + 1) * q_dim1 + 1], &c__1, &c_b507, &workd[*n + 1], &c__1, (ftnlen) 1);
	}

	/*     %----------------------------------------------------------% */
	/*     | Compute column 1 to kev of (V*Q) in backward order       | */
	/*     | taking advantage of the upper Hessenberg structure of Q. | */
	/*     %----------------------------------------------------------% */

	i__1 = *kev;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		i__2 = kplusp - i__ + 1;
		dgemv_ ("N", n, &i__2, &c_b348, &v[v_offset], ldv, &q[(*kev - i__ + 1) * q_dim1 + 1], &c__1, &c_b507, &workd[1], &c__1, (ftnlen) 1);
		dcopy_ (n, &workd[1], &c__1, &v[(kplusp - i__ + 1) * v_dim1 + 1], &c__1);
		/* L140: */
	}

	/*     %-------------------------------------------------% */
	/*     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). | */
	/*     %-------------------------------------------------% */

	dlacpy_ ("A", n, kev, &v[(kplusp - *kev + 1) * v_dim1 + 1], ldv, &v[v_offset], ldv, (ftnlen) 1);

	/*     %--------------------------------------------------------------% */
	/*     | Copy the (kev+1)-st column of (V*Q) in the appropriate place | */
	/*     %--------------------------------------------------------------% */

	if (h__[*kev + 1 + *kev * h_dim1] > 0.)
	{
		dcopy_ (n, &workd[*n + 1], &c__1, &v[(*kev + 1) * v_dim1 + 1], &c__1);
	}

	/*     %-------------------------------------% */
	/*     | Update the residual vector:         | */
	/*     |    r <- sigmak*r + betak*v(:,kev+1) | */
	/*     | where                               | */
	/*     |    sigmak = (e_{kplusp}'*Q)*e_{kev} | */
	/*     |    betak = e_{kev+1}'*H*e_{kev}     | */
	/*     %-------------------------------------% */

	dscal_ (n, &q[kplusp + *kev * q_dim1], &resid[1], &c__1);
	if (h__[*kev + 1 + *kev * h_dim1] > 0.)
	{
		daxpy_ (n, &h__[*kev + 1 + *kev * h_dim1], &v[(*kev + 1) * v_dim1 + 1], &c__1, &resid[1], &c__1);
	}

	if (msglvl > 1)
	{
		dvout_ (&debug_1.logfil, &c__1, &q[kplusp + *kev * q_dim1], &debug_1.ndigit, "_napps: sigmak = (e_{kev+p}^T*Q)*e_{kev}", (ftnlen) 40);
		dvout_ (&debug_1.logfil, &c__1, &h__[*kev + 1 + *kev * h_dim1], &debug_1.ndigit, "_napps: betak = e_{kev+1}^T*H*e_{kev}", (ftnlen) 37);
		ivout_ (&debug_1.logfil, &c__1, kev, &debug_1.ndigit, "_napps: Order " "of the final Hessenberg matrix ", (ftnlen) 45);
		if (msglvl > 2)
		{
			dmout_ (&debug_1.logfil, kev, kev, &h__[h_offset], ldh, &debug_1.ndigit, "_napps: updated Hessenberg matrix H for" " next iteration", (ftnlen) 54);
		}

	}

 L9000:
	second_ (&t1);
	timing_1.tnapps += t1 - t0;

	return 0;

	/*     %---------------% */
	/*     | End of dnapps | */
	/*     %---------------% */

}	/* dnapps_ */

/* \BeginDoc */

/* \Name: dnaup2 */

/* \Description: */
/*  Intermediate level interface called by dnaupd. */

/* \Usage: */
/*  call dnaup2 */
/*     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD, */
/*       ISHIFT, MXITER, V, LDV, H, LDH, RITZR, RITZI, BOUNDS, */
/*       Q, LDQ, WORKL, IPNTR, WORKD, INFO ) */

/* \Arguments */

/*  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in dnaupd. */
/*  MODE, ISHIFT, MXITER: see the definition of IPARAM in dnaupd. */

/*  NP      Integer.  (INPUT/OUTPUT) */
/*          Contains the number of implicit shifts to apply during */
/*          each Arnoldi iteration. */
/*          If ISHIFT=1, NP is adjusted dynamically at each iteration */
/*          to accelerate convergence and prevent stagnation. */
/*          This is also roughly equal to the number of matrix-vector */
/*          products (involving the operator OP) per Arnoldi iteration. */
/*          The logic for adjusting is contained within the current */
/*          subroutine. */
/*          If ISHIFT=0, NP is the number of shifts the user needs */
/*          to provide via reverse comunication. 0 < NP < NCV-NEV. */
/*          NP may be less than NCV-NEV for two reasons. The first, is */
/*          to keep complex conjugate pairs of "wanted" Ritz values */
/*          together. The second, is that a leading block of the current */
/*          upper Hessenberg matrix has split off and contains "unwanted" */
/*          Ritz values. */
/*          Upon termination of the IRA iteration, NP contains the number */
/*          of "converged" wanted Ritz values. */

/*  IUPD    Integer.  (INPUT) */
/*          IUPD .EQ. 0: use explicit restart instead implicit update. */
/*          IUPD .NE. 0: use implicit update. */

/*  V       Double precision N by (NEV+NP) array.  (INPUT/OUTPUT) */
/*          The Arnoldi basis vectors are returned in the first NEV */
/*          columns of V. */

/*  LDV     Integer.  (INPUT) */
/*          Leading dimension of V exactly as declared in the calling */
/*          program. */

/*  H       Double precision (NEV+NP) by (NEV+NP) array.  (OUTPUT) */
/*          H is used to store the generated upper Hessenberg matrix */

/*  LDH     Integer.  (INPUT) */
/*          Leading dimension of H exactly as declared in the calling */
/*          program. */

/*  RITZR,  Double precision arrays of length NEV+NP.  (OUTPUT) */
/*  RITZI   RITZR(1:NEV) (resp. RITZI(1:NEV)) contains the real (resp. */
/*          imaginary) part of the computed Ritz values of OP. */

/*  BOUNDS  Double precision array of length NEV+NP.  (OUTPUT) */
/*          BOUNDS(1:NEV) contain the error bounds corresponding to */
/*          the computed Ritz values. */

/*  Q       Double precision (NEV+NP) by (NEV+NP) array.  (WORKSPACE) */
/*          Private (replicated) work array used to accumulate the */
/*          rotation in the shift application step. */

/*  LDQ     Integer.  (INPUT) */
/*          Leading dimension of Q exactly as declared in the calling */
/*          program. */

/*  WORKL   Double precision work array of length at least */
/*          (NEV+NP)**2 + 3*(NEV+NP).  (INPUT/WORKSPACE) */
/*          Private (replicated) array on each PE or array allocated on */
/*          the front end.  It is used in shifts calculation, shifts */
/*          application and convergence checking. */

/*          On exit, the last 3*(NEV+NP) locations of WORKL contain */
/*          the Ritz values (real,imaginary) and associated Ritz */
/*          estimates of the current Hessenberg matrix.  They are */
/*          listed in the same order as returned from dneigh. */

/*          If ISHIFT .EQ. O and IDO .EQ. 3, the first 2*NP locations */
/*          of WORKL are used in reverse communication to hold the user */
/*          supplied shifts. */

/*  IPNTR   Integer array of length 3.  (OUTPUT) */
/*          Pointer to mark the starting locations in the WORKD for */
/*          vectors used by the Arnoldi iteration. */
/*          ------------------------------------------------------------- */
/*          IPNTR(1): pointer to the current operand vector X. */
/*          IPNTR(2): pointer to the current result vector Y. */
/*          IPNTR(3): pointer to the vector B * X when used in the */
/*                    shift-and-invert mode.  X is the current operand. */
/*          ------------------------------------------------------------- */

/*  WORKD   Double precision work array of length 3*N.  (WORKSPACE) */
/*          Distributed array to be used in the basic Arnoldi iteration */
/*          for reverse communication.  The user should not use WORKD */
/*          as temporary workspace during the iteration !!!!!!!!!! */
/*          See Data Distribution Note in DNAUPD. */

/*  INFO    Integer.  (INPUT/OUTPUT) */
/*          If INFO .EQ. 0, a randomly initial residual vector is used. */
/*          If INFO .NE. 0, RESID contains the initial residual vector, */
/*                          possibly from a previous run. */
/*          Error flag on output. */
/*          =     0: Normal return. */
/*          =     1: Maximum number of iterations taken. */
/*                   All possible eigenvalues of OP has been found. */
/*                   NP returns the number of converged Ritz values. */
/*          =     2: No shifts could be applied. */
/*          =    -8: Error return from LAPACK eigenvalue calculation; */
/*                   This should never happen. */
/*          =    -9: Starting vector is zero. */
/*          = -9999: Could not build an Arnoldi factorization. */
/*                   Size that was built in returned in NP. */

/* \EndDoc */

/* ----------------------------------------------------------------------- */

/* \BeginLib */

/* \Local variables: */
/*     xxxxxx  real */

/* \References: */
/*  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in */
/*     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992), */
/*     pp 357-385. */
/*  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly */
/*     Restarted Arnoldi Iteration", Rice University Technical Report */
/*     TR95-13, Department of Computational and Applied Mathematics. */

/* \Routines called: */
/*     dgetv0  ARPACK initial vector generation routine. */
/*     dnaitr  ARPACK Arnoldi factorization routine. */
/*     dnapps  ARPACK application of implicit shifts routine. */
/*     dnconv  ARPACK convergence of Ritz values routine. */
/*     dneigh  ARPACK compute Ritz values and error bounds routine. */
/*     dngets  ARPACK reorder Ritz values and error bounds routine. */
/*     dsortc  ARPACK sorting routine. */
/*     ivout   ARPACK utility routine that prints integers. */
/*     second  ARPACK utility routine for timing. */
/*     dmout   ARPACK utility routine that prints matrices */
/*     dvout   ARPACK utility routine that prints vectors. */
/*     dlamch  LAPACK routine that determines machine constants. */
/*     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully. */
/*     dcopy   Level 1 BLAS that copies one vector to another . */
/*     ddot    Level 1 BLAS that computes the scalar product of two vectors. */
/*     dnrm2   Level 1 BLAS that computes the norm of a vector. */
/*     dswap   Level 1 BLAS that swaps two vectors. */

/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \SCCS Information: @(#) */
/* FILE: naup2.F   SID: 2.8   DATE OF SID: 10/17/00   RELEASE: 2 */

/* \Remarks */
/*     1. None */

/* \EndLib */

/* ----------------------------------------------------------------------- */

static int
dnaup2_ (integer * ido, char *bmat, integer * n, char *which, integer * nev, integer * np, doublereal * tol, doublereal * resid,
			integer * mode, integer * iupd, integer * ishift, integer * mxiter,
			doublereal * v, integer * ldv, doublereal * h__, integer * ldh,
			doublereal * ritzr, doublereal * ritzi, doublereal * bounds, doublereal *
			q, integer * ldq, doublereal * workl, integer * ipntr, doublereal * workd, integer * info, ftnlen bmat_len, ftnlen which_len)
{
	/* System generated locals */
	integer h_dim1, h_offset, q_dim1, q_offset, v_dim1, v_offset, i__1, i__2;
	doublereal d__1, d__2;

	/* Builtin functions */
	double pow_dd (doublereal *, doublereal *);
	integer s_cmp (char *, char *, ftnlen, ftnlen);
	/* Subroutine */ int s_copy (char *, char *, ftnlen, ftnlen);
	double sqrt (doublereal);

	/* Local variables */
	static integer j;
	static real t0, t1, t2, t3;
	static integer kp[4], np0, nev0;
	static doublereal eps23;
	static integer ierr, iter;
	static doublereal temp;
	static logical getv0, cnorm;
	static integer nconv;
	static logical initv;
	static doublereal rnorm;
	static integer nevbef;
	static logical update;
	static char wprime[2];
	static logical ushift;
	static integer kplusp, msglvl, nptemp, numcnv;


	/*     %----------------------------------------------------% */
	/*     | Include files for debugging and timing information | */
	/*     %----------------------------------------------------% */


	/* \SCCS Information: @(#) */
	/* FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 */

	/*     %---------------------------------% */
	/*     | See debug.doc for documentation | */
	/*     %---------------------------------% */

	/*     %------------------% */
	/*     | Scalar Arguments | */
	/*     %------------------% */

	/*     %--------------------------------% */
	/*     | See stat.doc for documentation | */
	/*     %--------------------------------% */

	/* \SCCS Information: @(#) */
	/* FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 */



	/*     %-----------------% */
	/*     | Array Arguments | */
	/*     %-----------------% */


	/*     %------------% */
	/*     | Parameters | */
	/*     %------------% */


	/*     %---------------% */
	/*     | Local Scalars | */
	/*     %---------------% */


	/*     %-----------------------% */
	/*     | Local array arguments | */
	/*     %-----------------------% */


	/*     %----------------------% */
	/*     | External Subroutines | */
	/*     %----------------------% */


	/*     %--------------------% */
	/*     | External Functions | */
	/*     %--------------------% */


	/*     %---------------------% */
	/*     | Intrinsic Functions | */
	/*     %---------------------% */


	/*     %-----------------------% */
	/*     | Executable Statements | */
	/*     %-----------------------% */

	/* Parameter adjustments */
	--workd;
	--resid;
	--workl;
	--bounds;
	--ritzi;
	--ritzr;
	v_dim1 = *ldv;
	v_offset = 1 + v_dim1;
	v -= v_offset;
	h_dim1 = *ldh;
	h_offset = 1 + h_dim1;
	h__ -= h_offset;
	q_dim1 = *ldq;
	q_offset = 1 + q_dim1;
	q -= q_offset;
	--ipntr;

	/* Function Body */
	if (*ido == 0)
	{

		second_ (&t0);

		msglvl = debug_1.mnaup2;

		/*        %-------------------------------------% */
		/*        | Get the machine dependent constant. | */
		/*        %-------------------------------------% */

		eps23 = dlamch_ ("Epsilon-Machine", (ftnlen) 15);
		eps23 = pow_dd (&eps23, &c_b2616);

		nev0 = *nev;
		np0 = *np;

		/*        %-------------------------------------% */
		/*        | kplusp is the bound on the largest  | */
		/*        |        Lanczos factorization built. | */
		/*        | nconv is the current number of      | */
		/*        |        "converged" eigenvlues.      | */
		/*        | iter is the counter on the current  | */
		/*        |      iteration step.                | */
		/*        %-------------------------------------% */

		kplusp = *nev + *np;
		nconv = 0;
		iter = 0;

		/*        %---------------------------------------% */
		/*        | Set flags for computing the first NEV | */
		/*        | steps of the Arnoldi factorization.   | */
		/*        %---------------------------------------% */

		getv0 = TRUE_;
		update = FALSE_;
		ushift = FALSE_;
		cnorm = FALSE_;

		if (*info != 0)
		{

			/*           %--------------------------------------------% */
			/*           | User provides the initial residual vector. | */
			/*           %--------------------------------------------% */

			initv = TRUE_;
			*info = 0;
		}
		else
		{
			initv = FALSE_;
		}
	}

	/*     %---------------------------------------------% */
	/*     | Get a possibly random starting vector and   | */
	/*     | force it into the range of the operator OP. | */
	/*     %---------------------------------------------% */

	/* L10: */

	if (getv0)
	{
		dgetv0_ (ido, bmat, &c__1, &initv, n, &c__1, &v[v_offset], ldv, &resid[1], &rnorm, &ipntr[1], &workd[1], info, (ftnlen) 1);

		if (*ido != 99)
		{
			goto L9000;
		}

		if (rnorm == 0.)
		{

			/*           %-----------------------------------------% */
			/*           | The initial vector is zero. Error exit. | */
			/*           %-----------------------------------------% */

			*info = -9;
			goto L1100;
		}
		getv0 = FALSE_;
		*ido = 0;
	}

	/*     %-----------------------------------% */
	/*     | Back from reverse communication : | */
	/*     | continue with update step         | */
	/*     %-----------------------------------% */

	if (update)
	{
		goto L20;
	}

	/*     %-------------------------------------------% */
	/*     | Back from computing user specified shifts | */
	/*     %-------------------------------------------% */

	if (ushift)
	{
		goto L50;
	}

	/*     %-------------------------------------% */
	/*     | Back from computing residual norm   | */
	/*     | at the end of the current iteration | */
	/*     %-------------------------------------% */

	if (cnorm)
	{
		goto L100;
	}

	/*     %----------------------------------------------------------% */
	/*     | Compute the first NEV steps of the Arnoldi factorization | */
	/*     %----------------------------------------------------------% */

	dnaitr_ (ido, bmat, n, &c__0, nev, mode, &resid[1], &rnorm, &v[v_offset], ldv, &h__[h_offset], ldh, &ipntr[1], &workd[1], info, (ftnlen) 1);

	/*     %---------------------------------------------------% */
	/*     | ido .ne. 99 implies use of reverse communication  | */
	/*     | to compute operations involving OP and possibly B | */
	/*     %---------------------------------------------------% */

	if (*ido != 99)
	{
		goto L9000;
	}

	if (*info > 0)
	{
		*np = *info;
		*mxiter = iter;
		*info = -9999;
		goto L1200;
	}

	/*     %--------------------------------------------------------------% */
	/*     |                                                              | */
	/*     |           M A I N  ARNOLDI  I T E R A T I O N  L O O P       | */
	/*     |           Each iteration implicitly restarts the Arnoldi     | */
	/*     |           factorization in place.                            | */
	/*     |                                                              | */
	/*     %--------------------------------------------------------------% */

 L1000:

	++iter;

	if (msglvl > 0)
	{
		ivout_ (&debug_1.logfil, &c__1, &iter, &debug_1.ndigit, "_naup2: ****" " Start of major iteration number ****", (ftnlen) 49);
	}

	/*        %-----------------------------------------------------------% */
	/*        | Compute NP additional steps of the Arnoldi factorization. | */
	/*        | Adjust NP since NEV might have been updated by last call  | */
	/*        | to the shift application routine dnapps.                  | */
	/*        %-----------------------------------------------------------% */

	*np = kplusp - *nev;

	if (msglvl > 1)
	{
		ivout_ (&debug_1.logfil, &c__1, nev, &debug_1.ndigit, "_naup2: The le" "ngth of the current Arnoldi factorization", (ftnlen) 55);
		ivout_ (&debug_1.logfil, &c__1, np, &debug_1.ndigit, "_naup2: Extend " "the Arnoldi factorization by", (ftnlen) 43);
	}

	/*        %-----------------------------------------------------------% */
	/*        | Compute NP additional steps of the Arnoldi factorization. | */
	/*        %-----------------------------------------------------------% */

	*ido = 0;
 L20:
	update = TRUE_;

	dnaitr_ (ido, bmat, n, nev, np, mode, &resid[1], &rnorm, &v[v_offset], ldv, &h__[h_offset], ldh, &ipntr[1], &workd[1], info, (ftnlen) 1);

	/*        %---------------------------------------------------% */
	/*        | ido .ne. 99 implies use of reverse communication  | */
	/*        | to compute operations involving OP and possibly B | */
	/*        %---------------------------------------------------% */

	if (*ido != 99)
	{
		goto L9000;
	}

	if (*info > 0)
	{
		*np = *info;
		*mxiter = iter;
		*info = -9999;
		goto L1200;
	}
	update = FALSE_;

	if (msglvl > 1)
	{
		dvout_ (&debug_1.logfil, &c__1, &rnorm, &debug_1.ndigit, "_naup2: Cor" "responding B-norm of the residual", (ftnlen) 44);
	}

	/*        %--------------------------------------------------------% */
	/*        | Compute the eigenvalues and corresponding error bounds | */
	/*        | of the current upper Hessenberg matrix.                | */
	/*        %--------------------------------------------------------% */

	dneigh_ (&rnorm, &kplusp, &h__[h_offset], ldh, &ritzr[1], &ritzi[1], &bounds[1], &q[q_offset], ldq, &workl[1], &ierr);

	if (ierr != 0)
	{
		*info = -8;
		goto L1200;
	}

	/*        %----------------------------------------------------% */
	/*        | Make a copy of eigenvalues and corresponding error | */
	/*        | bounds obtained from dneigh.                       | */
	/*        %----------------------------------------------------% */

	/* Computing 2nd power */
	i__1 = kplusp;
	dcopy_ (&kplusp, &ritzr[1], &c__1, &workl[i__1 * i__1 + 1], &c__1);
	/* Computing 2nd power */
	i__1 = kplusp;
	dcopy_ (&kplusp, &ritzi[1], &c__1, &workl[i__1 * i__1 + kplusp + 1], &c__1);
	/* Computing 2nd power */
	i__1 = kplusp;
	dcopy_ (&kplusp, &bounds[1], &c__1, &workl[i__1 * i__1 + (kplusp << 1) + 1], &c__1);

	/*        %---------------------------------------------------% */
	/*        | Select the wanted Ritz values and their bounds    | */
	/*        | to be used in the convergence test.               | */
	/*        | The wanted part of the spectrum and corresponding | */
	/*        | error bounds are in the last NEV loc. of RITZR,   | */
	/*        | RITZI and BOUNDS respectively. The variables NEV  | */
	/*        | and NP may be updated if the NEV-th wanted Ritz   | */
	/*        | value has a non zero imaginary part. In this case | */
	/*        | NEV is increased by one and NP decreased by one.  | */
	/*        | NOTE: The last two arguments of dngets are no     | */
	/*        | longer used as of version 2.1.                    | */
	/*        %---------------------------------------------------% */

	*nev = nev0;
	*np = np0;
	numcnv = *nev;
	dngets_ (ishift, which, nev, np, &ritzr[1], &ritzi[1], &bounds[1], &workl[1], &workl[*np + 1], (ftnlen) 2);
	if (*nev == nev0 + 1)
	{
		numcnv = nev0 + 1;
	}

	/*        %-------------------% */
	/*        | Convergence test. | */
	/*        %-------------------% */

	dcopy_ (nev, &bounds[*np + 1], &c__1, &workl[(*np << 1) + 1], &c__1);
	dnconv_ (nev, &ritzr[*np + 1], &ritzi[*np + 1], &workl[(*np << 1) + 1], tol, &nconv);

	if (msglvl > 2)
	{
		kp[0] = *nev;
		kp[1] = *np;
		kp[2] = numcnv;
		kp[3] = nconv;
		ivout_ (&debug_1.logfil, &c__4, kp, &debug_1.ndigit, "_naup2: NEV, NP" ", NUMCNV, NCONV are", (ftnlen) 34);
		dvout_ (&debug_1.logfil, &kplusp, &ritzr[1], &debug_1.ndigit, "_naup2" ": Real part of the eigenvalues of H", (ftnlen) 41);
		dvout_ (&debug_1.logfil, &kplusp, &ritzi[1], &debug_1.ndigit, "_naup2" ": Imaginary part of the eigenvalues of H", (ftnlen) 46);
		dvout_ (&debug_1.logfil, &kplusp, &bounds[1], &debug_1.ndigit, "_naup" "2: Ritz estimates of the current NCV Ritz values", (ftnlen) 53);
	}

	/*        %---------------------------------------------------------% */
	/*        | Count the number of unwanted Ritz values that have zero | */
	/*        | Ritz estimates. If any Ritz estimates are equal to zero | */
	/*        | then a leading block of H of order equal to at least    | */
	/*        | the number of Ritz values with zero Ritz estimates has  | */
	/*        | split off. None of these Ritz values may be removed by  | */
	/*        | shifting. Decrease NP the number of shifts to apply. If | */
	/*        | no shifts may be applied, then prepare to exit          | */
	/*        %---------------------------------------------------------% */

	nptemp = *np;
	i__1 = nptemp;
	for (j = 1; j <= i__1; ++j)
	{
		if (bounds[j] == 0.)
		{
			--(*np);
			++(*nev);
		}
		/* L30: */
	}

	if (nconv >= numcnv || iter > *mxiter || *np == 0)
	{

		if (msglvl > 4)
		{
			/* Computing 2nd power */
			i__1 = kplusp;
			dvout_ (&debug_1.logfil, &kplusp, &workl[i__1 * i__1 + 1], &debug_1.ndigit, "_naup2: Real part of the eig computed b" "y _neigh:", (ftnlen) 48);
			/* Computing 2nd power */
			i__1 = kplusp;
			dvout_ (&debug_1.logfil, &kplusp, &workl[i__1 * i__1 + kplusp + 1],
					  &debug_1.ndigit, "_naup2: Imag part of the eig computed" " by _neigh:", (ftnlen) 48);
			/* Computing 2nd power */
			i__1 = kplusp;
			dvout_ (&debug_1.logfil, &kplusp, &workl[i__1 * i__1 + (kplusp <<
																					  1) + 1], &debug_1.ndigit, "_naup2: Ritz eistmates comput" "ed by _neigh:", (ftnlen) 42);
		}

		/*           %------------------------------------------------% */
		/*           | Prepare to exit. Put the converged Ritz values | */
		/*           | and corresponding bounds in RITZ(1:NCONV) and  | */
		/*           | BOUNDS(1:NCONV) respectively. Then sort. Be    | */
		/*           | careful when NCONV > NP                        | */
		/*           %------------------------------------------------% */

		/*           %------------------------------------------% */
		/*           |  Use h( 3,1 ) as storage to communicate  | */
		/*           |  rnorm to _neupd if needed               | */
		/*           %------------------------------------------% */
		h__[h_dim1 + 3] = rnorm;

		/*           %----------------------------------------------% */
		/*           | To be consistent with dngets, we first do a  | */
		/*           | pre-processing sort in order to keep complex | */
		/*           | conjugate pairs together.  This is similar   | */
		/*           | to the pre-processing sort used in dngets    | */
		/*           | except that the sort is done in the opposite | */
		/*           | order.                                       | */
		/*           %----------------------------------------------% */

		if (s_cmp (which, "LM", (ftnlen) 2, (ftnlen) 2) == 0)
		{
			s_copy (wprime, "SR", (ftnlen) 2, (ftnlen) 2);
		}
		if (s_cmp (which, "SM", (ftnlen) 2, (ftnlen) 2) == 0)
		{
			s_copy (wprime, "LR", (ftnlen) 2, (ftnlen) 2);
		}
		if (s_cmp (which, "LR", (ftnlen) 2, (ftnlen) 2) == 0)
		{
			s_copy (wprime, "SM", (ftnlen) 2, (ftnlen) 2);
		}
		if (s_cmp (which, "SR", (ftnlen) 2, (ftnlen) 2) == 0)
		{
			s_copy (wprime, "LM", (ftnlen) 2, (ftnlen) 2);
		}
		if (s_cmp (which, "LI", (ftnlen) 2, (ftnlen) 2) == 0)
		{
			s_copy (wprime, "SM", (ftnlen) 2, (ftnlen) 2);
		}
		if (s_cmp (which, "SI", (ftnlen) 2, (ftnlen) 2) == 0)
		{
			s_copy (wprime, "LM", (ftnlen) 2, (ftnlen) 2);
		}

		dsortc_ (wprime, &c_true, &kplusp, &ritzr[1], &ritzi[1], &bounds[1], (ftnlen) 2);

		/*           %----------------------------------------------% */
		/*           | Now sort Ritz values so that converged Ritz  | */
		/*           | values appear within the first NEV locations | */
		/*           | of ritzr, ritzi and bounds, and the most     | */
		/*           | desired one appears at the front.            | */
		/*           %----------------------------------------------% */

		if (s_cmp (which, "LM", (ftnlen) 2, (ftnlen) 2) == 0)
		{
			s_copy (wprime, "SM", (ftnlen) 2, (ftnlen) 2);
		}
		if (s_cmp (which, "SM", (ftnlen) 2, (ftnlen) 2) == 0)
		{
			s_copy (wprime, "LM", (ftnlen) 2, (ftnlen) 2);
		}
		if (s_cmp (which, "LR", (ftnlen) 2, (ftnlen) 2) == 0)
		{
			s_copy (wprime, "SR", (ftnlen) 2, (ftnlen) 2);
		}
		if (s_cmp (which, "SR", (ftnlen) 2, (ftnlen) 2) == 0)
		{
			s_copy (wprime, "LR", (ftnlen) 2, (ftnlen) 2);
		}
		if (s_cmp (which, "LI", (ftnlen) 2, (ftnlen) 2) == 0)
		{
			s_copy (wprime, "SI", (ftnlen) 2, (ftnlen) 2);
		}
		if (s_cmp (which, "SI", (ftnlen) 2, (ftnlen) 2) == 0)
		{
			s_copy (wprime, "LI", (ftnlen) 2, (ftnlen) 2);
		}

		dsortc_ (wprime, &c_true, &kplusp, &ritzr[1], &ritzi[1], &bounds[1], (ftnlen) 2);

		/*           %--------------------------------------------------% */
		/*           | Scale the Ritz estimate of each Ritz value       | */
		/*           | by 1 / max(eps23,magnitude of the Ritz value).   | */
		/*           %--------------------------------------------------% */

		i__1 = numcnv;
		for (j = 1; j <= i__1; ++j)
		{
			/* Computing MAX */
			d__1 = eps23, d__2 = dlapy2_ (&ritzr[j], &ritzi[j]);
			temp = max (d__1, d__2);
			bounds[j] /= temp;
			/* L35: */
		}

		/*           %----------------------------------------------------% */
		/*           | Sort the Ritz values according to the scaled Ritz  | */
		/*           | esitmates.  This will push all the converged ones  | */
		/*           | towards the front of ritzr, ritzi, bounds          | */
		/*           | (in the case when NCONV < NEV.)                    | */
		/*           %----------------------------------------------------% */

		s_copy (wprime, "LR", (ftnlen) 2, (ftnlen) 2);
		dsortc_ (wprime, &c_true, &numcnv, &bounds[1], &ritzr[1], &ritzi[1], (ftnlen) 2);

		/*           %----------------------------------------------% */
		/*           | Scale the Ritz estimate back to its original | */
		/*           | value.                                       | */
		/*           %----------------------------------------------% */

		i__1 = numcnv;
		for (j = 1; j <= i__1; ++j)
		{
			/* Computing MAX */
			d__1 = eps23, d__2 = dlapy2_ (&ritzr[j], &ritzi[j]);
			temp = max (d__1, d__2);
			bounds[j] *= temp;
			/* L40: */
		}

		/*           %------------------------------------------------% */
		/*           | Sort the converged Ritz values again so that   | */
		/*           | the "threshold" value appears at the front of  | */
		/*           | ritzr, ritzi and bound.                        | */
		/*           %------------------------------------------------% */

		dsortc_ (which, &c_true, &nconv, &ritzr[1], &ritzi[1], &bounds[1], (ftnlen) 2);

		if (msglvl > 1)
		{
			dvout_ (&debug_1.logfil, &kplusp, &ritzr[1], &debug_1.ndigit, "_naup2: Sorted real part of the eigenvalues", (ftnlen) 43);
			dvout_ (&debug_1.logfil, &kplusp, &ritzi[1], &debug_1.ndigit, "_naup2: Sorted imaginary part of the eigenvalues", (ftnlen) 48);
			dvout_ (&debug_1.logfil, &kplusp, &bounds[1], &debug_1.ndigit, "_naup2: Sorted ritz estimates.", (ftnlen) 30);
		}

		/*           %------------------------------------% */
		/*           | Max iterations have been exceeded. | */
		/*           %------------------------------------% */

		if (iter > *mxiter && nconv < numcnv)
		{
			*info = 1;
		}

		/*           %---------------------% */
		/*           | No shifts to apply. | */
		/*           %---------------------% */

		if (*np == 0 && nconv < numcnv)
		{
			*info = 2;
		}

		*np = nconv;
		goto L1100;

	}
	else if (nconv < numcnv && *ishift == 1)
	{

		/*           %-------------------------------------------------% */
		/*           | Do not have all the requested eigenvalues yet.  | */
		/*           | To prevent possible stagnation, adjust the size | */
		/*           | of NEV.                                         | */
		/*           %-------------------------------------------------% */

		nevbef = *nev;
		/* Computing MIN */
		i__1 = nconv, i__2 = *np / 2;
		*nev += min (i__1, i__2);
		if (*nev == 1 && kplusp >= 6)
		{
			*nev = kplusp / 2;
		}
		else if (*nev == 1 && kplusp > 3)
		{
			*nev = 2;
		}
		*np = kplusp - *nev;

		/*           %---------------------------------------% */
		/*           | If the size of NEV was just increased | */
		/*           | resort the eigenvalues.               | */
		/*           %---------------------------------------% */

		if (nevbef < *nev)
		{
			dngets_ (ishift, which, nev, np, &ritzr[1], &ritzi[1], &bounds[1], &workl[1], &workl[*np + 1], (ftnlen) 2);
		}

	}

	if (msglvl > 0)
	{
		ivout_ (&debug_1.logfil, &c__1, &nconv, &debug_1.ndigit, "_naup2: no." " of \"converged\" Ritz values at this iter.", (ftnlen) 52);
		if (msglvl > 1)
		{
			kp[0] = *nev;
			kp[1] = *np;
			ivout_ (&debug_1.logfil, &c__2, kp, &debug_1.ndigit, "_naup2: NEV" " and NP are", (ftnlen) 22);
			dvout_ (&debug_1.logfil, nev, &ritzr[*np + 1], &debug_1.ndigit, "_naup2: \"wanted\" Ritz values -- real part", (ftnlen) 41);
			dvout_ (&debug_1.logfil, nev, &ritzi[*np + 1], &debug_1.ndigit, "_naup2: \"wanted\" Ritz values -- imag part", (ftnlen) 41);
			dvout_ (&debug_1.logfil, nev, &bounds[*np + 1], &debug_1.ndigit, "_naup2: Ritz estimates of the \"wanted\" values ", (ftnlen) 46);
		}
	}

	if (*ishift == 0)
	{

		/*           %-------------------------------------------------------% */
		/*           | User specified shifts: reverse comminucation to       | */
		/*           | compute the shifts. They are returned in the first    | */
		/*           | 2*NP locations of WORKL.                              | */
		/*           %-------------------------------------------------------% */

		ushift = TRUE_;
		*ido = 3;
		goto L9000;
	}

 L50:

	/*        %------------------------------------% */
	/*        | Back from reverse communication;   | */
	/*        | User specified shifts are returned | */
	/*        | in WORKL(1:2*NP)                   | */
	/*        %------------------------------------% */

	ushift = FALSE_;

	if (*ishift == 0)
	{

		/*            %----------------------------------% */
		/*            | Move the NP shifts from WORKL to | */
		/*            | RITZR, RITZI to free up WORKL    | */
		/*            | for non-exact shift case.        | */
		/*            %----------------------------------% */

		dcopy_ (np, &workl[1], &c__1, &ritzr[1], &c__1);
		dcopy_ (np, &workl[*np + 1], &c__1, &ritzi[1], &c__1);
	}

	if (msglvl > 2)
	{
		ivout_ (&debug_1.logfil, &c__1, np, &debug_1.ndigit, "_naup2: The num" "ber of shifts to apply ", (ftnlen) 38);
		dvout_ (&debug_1.logfil, np, &ritzr[1], &debug_1.ndigit, "_naup2: Rea" "l part of the shifts", (ftnlen) 31);
		dvout_ (&debug_1.logfil, np, &ritzi[1], &debug_1.ndigit, "_naup2: Ima" "ginary part of the shifts", (ftnlen) 36);
		if (*ishift == 1)
		{
			dvout_ (&debug_1.logfil, np, &bounds[1], &debug_1.ndigit, "_naup2" ": Ritz estimates of the shifts", (ftnlen) 36);
		}
	}

	/*        %---------------------------------------------------------% */
	/*        | Apply the NP implicit shifts by QR bulge chasing.       | */
	/*        | Each shift is applied to the whole upper Hessenberg     | */
	/*        | matrix H.                                               | */
	/*        | The first 2*N locations of WORKD are used as workspace. | */
	/*        %---------------------------------------------------------% */

	dnapps_ (n, nev, np, &ritzr[1], &ritzi[1], &v[v_offset], ldv, &h__[h_offset], ldh, &resid[1], &q[q_offset], ldq, &workl[1], &workd[1]);

	/*        %---------------------------------------------% */
	/*        | Compute the B-norm of the updated residual. | */
	/*        | Keep B*RESID in WORKD(1:N) to be used in    | */
	/*        | the first step of the next call to dnaitr.  | */
	/*        %---------------------------------------------% */

	cnorm = TRUE_;
	second_ (&t2);
	if (*(unsigned char *) bmat == 'G')
	{
		++timing_1.nbx;
		dcopy_ (n, &resid[1], &c__1, &workd[*n + 1], &c__1);
		ipntr[1] = *n + 1;
		ipntr[2] = 1;
		*ido = 2;

		/*           %----------------------------------% */
		/*           | Exit in order to compute B*RESID | */
		/*           %----------------------------------% */

		goto L9000;
	}
	else if (*(unsigned char *) bmat == 'I')
	{
		dcopy_ (n, &resid[1], &c__1, &workd[1], &c__1);
	}

 L100:

	/*        %----------------------------------% */
	/*        | Back from reverse communication; | */
	/*        | WORKD(1:N) := B*RESID            | */
	/*        %----------------------------------% */

	if (*(unsigned char *) bmat == 'G')
	{
		second_ (&t3);
		timing_1.tmvbx += t3 - t2;
	}

	if (*(unsigned char *) bmat == 'G')
	{
		rnorm = ddot_ (n, &resid[1], &c__1, &workd[1], &c__1);
		rnorm = sqrt ((abs (rnorm)));
	}
	else if (*(unsigned char *) bmat == 'I')
	{
		rnorm = dnrm2_ (n, &resid[1], &c__1);
	}
	cnorm = FALSE_;

	if (msglvl > 2)
	{
		dvout_ (&debug_1.logfil, &c__1, &rnorm, &debug_1.ndigit, "_naup2: B-n" "orm of residual for compressed factorization", (ftnlen) 55);
		dmout_ (&debug_1.logfil, nev, nev, &h__[h_offset], ldh, &debug_1.ndigit, "_naup2: Compressed upper Hessenberg matrix H", (ftnlen) 44);
	}

	goto L1000;

	/*     %---------------------------------------------------------------% */
	/*     |                                                               | */
	/*     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  | */
	/*     |                                                               | */
	/*     %---------------------------------------------------------------% */

 L1100:

	*mxiter = iter;
	*nev = numcnv;

 L1200:
	*ido = 99;

	/*     %------------% */
	/*     | Error Exit | */
	/*     %------------% */

	second_ (&t1);
	timing_1.tnaup2 = t1 - t0;

 L9000:

	/*     %---------------% */
	/*     | End of dnaup2 | */
	/*     %---------------% */

	return 0;
}	/* dnaup2_ */

/* \BeginDoc */

/* \Name: dnaupd */

/* \Description: */
/*  Reverse communication interface for the Implicitly Restarted Arnoldi */
/*  iteration. This subroutine computes approximations to a few eigenpairs */
/*  of a linear operator "OP" with respect to a semi-inner product defined by */
/*  a symmetric positive semi-definite real matrix B. B may be the identity */
/*  matrix. NOTE: If the linear operator "OP" is real and symmetric */
/*  with respect to the real positive semi-definite symmetric matrix B, */
/*  i.e. B*OP = (OP`)*B, then subroutine dsaupd should be used instead. */

/*  The computed approximate eigenvalues are called Ritz values and */
/*  the corresponding approximate eigenvectors are called Ritz vectors. */

/*  dnaupd is usually called iteratively to solve one of the */
/*  following problems: */

/*  Mode 1:  A*x = lambda*x. */
/*           ===> OP = A  and  B = I. */

/*  Mode 2:  A*x = lambda*M*x, M symmetric positive definite */
/*           ===> OP = inv[M]*A  and  B = M. */
/*           ===> (If M can be factored see remark 3 below) */

/*  Mode 3:  A*x = lambda*M*x, M symmetric semi-definite */
/*           ===> OP = Real_Part{ inv[A - sigma*M]*M }  and  B = M. */
/*           ===> shift-and-invert mode (in real arithmetic) */
/*           If OP*x = amu*x, then */
/*           amu = 1/2 * [ 1/(lambda-sigma) + 1/(lambda-conjg(sigma)) ]. */
/*           Note: If sigma is real, i.e. imaginary part of sigma is zero; */
/*                 Real_Part{ inv[A - sigma*M]*M } == inv[A - sigma*M]*M */
/*                 amu == 1/(lambda-sigma). */

/*  Mode 4:  A*x = lambda*M*x, M symmetric semi-definite */
/*           ===> OP = Imaginary_Part{ inv[A - sigma*M]*M }  and  B = M. */
/*           ===> shift-and-invert mode (in real arithmetic) */
/*           If OP*x = amu*x, then */
/*           amu = 1/2i * [ 1/(lambda-sigma) - 1/(lambda-conjg(sigma)) ]. */

/*  Both mode 3 and 4 give the same enhancement to eigenvalues close to */
/*  the (complex) shift sigma.  However, as lambda goes to infinity, */
/*  the operator OP in mode 4 dampens the eigenvalues more strongly than */
/*  does OP defined in mode 3. */

/*  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v */
/*        should be accomplished either by a direct method */
/*        using a sparse matrix factorization and solving */

/*           [A - sigma*M]*w = v  or M*w = v, */

/*        or through an iterative method for solving these */
/*        systems.  If an iterative method is used, the */
/*        convergence test must be more stringent than */
/*        the accuracy requirements for the eigenvalue */
/*        approximations. */

/* \Usage: */
/*  call dnaupd */
/*     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, */
/*       IPNTR, WORKD, WORKL, LWORKL, INFO ) */

/* \Arguments */
/*  IDO     Integer.  (INPUT/OUTPUT) */
/*          Reverse communication flag.  IDO must be zero on the first */
/*          call to dnaupd.  IDO will be set internally to */
/*          indicate the type of operation to be performed.  Control is */
/*          then given back to the calling routine which has the */
/*          responsibility to carry out the requested operation and call */
/*          dnaupd with the result.  The operand is given in */
/*          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)). */
/*          ------------------------------------------------------------- */
/*          IDO =  0: first call to the reverse communication interface */
/*          IDO = -1: compute  Y = OP * X  where */
/*                    IPNTR(1) is the pointer into WORKD for X, */
/*                    IPNTR(2) is the pointer into WORKD for Y. */
/*                    This is for the initialization phase to force the */
/*                    starting vector into the range of OP. */
/*          IDO =  1: compute  Y = OP * X  where */
/*                    IPNTR(1) is the pointer into WORKD for X, */
/*                    IPNTR(2) is the pointer into WORKD for Y. */
/*                    In mode 3 and 4, the vector B * X is already */
/*                    available in WORKD(ipntr(3)).  It does not */
/*                    need to be recomputed in forming OP * X. */
/*          IDO =  2: compute  Y = B * X  where */
/*                    IPNTR(1) is the pointer into WORKD for X, */
/*                    IPNTR(2) is the pointer into WORKD for Y. */
/*          IDO =  3: compute the IPARAM(8) real and imaginary parts */
/*                    of the shifts where INPTR(14) is the pointer */
/*                    into WORKL for placing the shifts. See Remark */
/*                    5 below. */
/*          IDO = 99: done */
/*          ------------------------------------------------------------- */

/*  BMAT    Character*1.  (INPUT) */
/*          BMAT specifies the type of the matrix B that defines the */
/*          semi-inner product for the operator OP. */
/*          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x */
/*          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*B*x */

/*  N       Integer.  (INPUT) */
/*          Dimension of the eigenproblem. */

/*  WHICH   Character*2.  (INPUT) */
/*          'LM' -> want the NEV eigenvalues of largest magnitude. */
/*          'SM' -> want the NEV eigenvalues of smallest magnitude. */
/*          'LR' -> want the NEV eigenvalues of largest real part. */
/*          'SR' -> want the NEV eigenvalues of smallest real part. */
/*          'LI' -> want the NEV eigenvalues of largest imaginary part. */
/*          'SI' -> want the NEV eigenvalues of smallest imaginary part. */

/*  NEV     Integer.  (INPUT/OUTPUT) */
/*          Number of eigenvalues of OP to be computed. 0 < NEV < N-1. */

/*  TOL     Double precision scalar.  (INPUT) */
/*          Stopping criterion: the relative accuracy of the Ritz value */
/*          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)) */
/*          where ABS(RITZ(I)) is the magnitude when RITZ(I) is complex. */
/*          DEFAULT = DLAMCH('EPS')  (machine precision as computed */
/*                    by the LAPACK auxiliary subroutine DLAMCH). */

/*  RESID   Double precision array of length N.  (INPUT/OUTPUT) */
/*          On INPUT: */
/*          If INFO .EQ. 0, a random initial residual vector is used. */
/*          If INFO .NE. 0, RESID contains the initial residual vector, */
/*                          possibly from a previous run. */
/*          On OUTPUT: */
/*          RESID contains the final residual vector. */

/*  NCV     Integer.  (INPUT) */
/*          Number of columns of the matrix V. NCV must satisfy the two */
/*          inequalities 2 <= NCV-NEV and NCV <= N. */
/*          This will indicate how many Arnoldi vectors are generated */
/*          at each iteration.  After the startup phase in which NEV */
/*          Arnoldi vectors are generated, the algorithm generates */
/*          approximately NCV-NEV Arnoldi vectors at each subsequent update */
/*          iteration. Most of the cost in generating each Arnoldi vector is */
/*          in the matrix-vector operation OP*x. */
/*          NOTE: 2 <= NCV-NEV in order that complex conjugate pairs of Ritz */
/*          values are kept together. (See remark 4 below) */

/*  V       Double precision array N by NCV.  (OUTPUT) */
/*          Contains the final set of Arnoldi basis vectors. */

/*  LDV     Integer.  (INPUT) */
/*          Leading dimension of V exactly as declared in the calling program. */

/*  IPARAM  Integer array of length 11.  (INPUT/OUTPUT) */
/*          IPARAM(1) = ISHIFT: method for selecting the implicit shifts. */
/*          The shifts selected at each iteration are used to restart */
/*          the Arnoldi iteration in an implicit fashion. */
/*          ------------------------------------------------------------- */
/*          ISHIFT = 0: the shifts are provided by the user via */
/*                      reverse communication.  The real and imaginary */
/*                      parts of the NCV eigenvalues of the Hessenberg */
/*                      matrix H are returned in the part of the WORKL */
/*                      array corresponding to RITZR and RITZI. See remark */
/*                      5 below. */
/*          ISHIFT = 1: exact shifts with respect to the current */
/*                      Hessenberg matrix H.  This is equivalent to */
/*                      restarting the iteration with a starting vector */
/*                      that is a linear combination of approximate Schur */
/*                      vectors associated with the "wanted" Ritz values. */
/*          ------------------------------------------------------------- */

/*          IPARAM(2) = No longer referenced. */

/*          IPARAM(3) = MXITER */
/*          On INPUT:  maximum number of Arnoldi update iterations allowed. */
/*          On OUTPUT: actual number of Arnoldi update iterations taken. */

/*          IPARAM(4) = NB: blocksize to be used in the recurrence. */
/*          The code currently works only for NB = 1. */

/*          IPARAM(5) = NCONV: number of "converged" Ritz values. */
/*          This represents the number of Ritz values that satisfy */
/*          the convergence criterion. */

/*          IPARAM(6) = IUPD */
/*          No longer referenced. Implicit restarting is ALWAYS used. */

/*          IPARAM(7) = MODE */
/*          On INPUT determines what type of eigenproblem is being solved. */
/*          Must be 1,2,3,4; See under \Description of dnaupd for the */
/*          four modes available. */

/*          IPARAM(8) = NP */
/*          When ido = 3 and the user provides shifts through reverse */
/*          communication (IPARAM(1)=0), dnaupd returns NP, the number */
/*          of shifts the user is to provide. 0 < NP <=NCV-NEV. See Remark */
/*          5 below. */

/*          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO, */
/*          OUTPUT: NUMOP  = total number of OP*x operations, */
/*                  NUMOPB = total number of B*x operations if BMAT='G', */
/*                  NUMREO = total number of steps of re-orthogonalization. */

/*  IPNTR   Integer array of length 14.  (OUTPUT) */
/*          Pointer to mark the starting locations in the WORKD and WORKL */
/*          arrays for matrices/vectors used by the Arnoldi iteration. */
/*          ------------------------------------------------------------- */
/*          IPNTR(1): pointer to the current operand vector X in WORKD. */
/*          IPNTR(2): pointer to the current result vector Y in WORKD. */
/*          IPNTR(3): pointer to the vector B * X in WORKD when used in */
/*                    the shift-and-invert mode. */
/*          IPNTR(4): pointer to the next available location in WORKL */
/*                    that is untouched by the program. */
/*          IPNTR(5): pointer to the NCV by NCV upper Hessenberg matrix */
/*                    H in WORKL. */
/*          IPNTR(6): pointer to the real part of the ritz value array */
/*                    RITZR in WORKL. */
/*          IPNTR(7): pointer to the imaginary part of the ritz value array */
/*                    RITZI in WORKL. */
/*          IPNTR(8): pointer to the Ritz estimates in array WORKL associated */
/*                    with the Ritz values located in RITZR and RITZI in WORKL. */

/*          IPNTR(14): pointer to the NP shifts in WORKL. See Remark 5 below. */

/*          Note: IPNTR(9:13) is only referenced by dneupd. See Remark 2 below. */

/*          IPNTR(9):  pointer to the real part of the NCV RITZ values of the */
/*                     original system. */
/*          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of */
/*                     the original system. */
/*          IPNTR(11): pointer to the NCV corresponding error bounds. */
/*          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular */
/*                     Schur matrix for H. */
/*          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors */
/*                     of the upper Hessenberg matrix H. Only referenced by */
/*                     dneupd if RVEC = .TRUE. See Remark 2 below. */
/*          ------------------------------------------------------------- */

/*  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION) */
/*          Distributed array to be used in the basic Arnoldi iteration */
/*          for reverse communication.  The user should not use WORKD */
/*          as temporary workspace during the iteration. Upon termination */
/*          WORKD(1:N) contains B*RESID(1:N). If an invariant subspace */
/*          associated with the converged Ritz values is desired, see remark */
/*          2 below, subroutine dneupd uses this output. */
/*          See Data Distribution Note below. */

/*  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE) */
/*          Private (replicated) array on each PE or array allocated on */
/*          the front end.  See Data Distribution Note below. */

/*  LWORKL  Integer.  (INPUT) */
/*          LWORKL must be at least 3*NCV**2 + 6*NCV. */

/*  INFO    Integer.  (INPUT/OUTPUT) */
/*          If INFO .EQ. 0, a randomly initial residual vector is used. */
/*          If INFO .NE. 0, RESID contains the initial residual vector, */
/*                          possibly from a previous run. */
/*          Error flag on output. */
/*          =  0: Normal exit. */
/*          =  1: Maximum number of iterations taken. */
/*                All possible eigenvalues of OP has been found. IPARAM(5) */
/*                returns the number of wanted converged Ritz values. */
/*          =  2: No longer an informational error. Deprecated starting */
/*                with release 2 of ARPACK. */
/*          =  3: No shifts could be applied during a cycle of the */
/*                Implicitly restarted Arnoldi iteration. One possibility */
/*                is to increase the size of NCV relative to NEV. */
/*                See remark 4 below. */
/*          = -1: N must be positive. */
/*          = -2: NEV must be positive. */
/*          = -3: NCV-NEV >= 2 and less than or equal to N. */
/*          = -4: The maximum number of Arnoldi update iteration */
/*                must be greater than zero. */
/*          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI' */
/*          = -6: BMAT must be one of 'I' or 'G'. */
/*          = -7: Length of private work array is not sufficient. */
/*          = -8: Error return from LAPACK eigenvalue calculation; */
/*          = -9: Starting vector is zero. */
/*          = -10: IPARAM(7) must be 1,2,3,4. */
/*          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable. */
/*          = -12: IPARAM(1) must be equal to 0 or 1. */
/*          = -9999: Could not build an Arnoldi factorization. */
/*                   IPARAM(5) returns the size of the current Arnoldi */
/*                   factorization. */

/* \Remarks */
/*  1. The computed Ritz values are approximate eigenvalues of OP. The */
/*     selection of WHICH should be made with this in mind when */
/*     Mode = 3 and 4.  After convergence, approximate eigenvalues of the */
/*     original problem may be obtained with the ARPACK subroutine dneupd. */

/*  2. If a basis for the invariant subspace corresponding to the converged Ritz */
/*     values is needed, the user must call dneupd immediately following */
/*     completion of dnaupd. This is new starting with release 2 of ARPACK. */

/*  3. If M can be factored into a Cholesky factorization M = LL` */
/*     then Mode = 2 should not be selected.  Instead one should use */
/*     Mode = 1 with  OP = inv(L)*A*inv(L`).  Appropriate triangular */
/*     linear systems should be solved with L and L` rather */
/*     than computing inverses.  After convergence, an approximate */
/*     eigenvector z of the original problem is recovered by solving */
/*     L`z = x  where x is a Ritz vector of OP. */

/*  4. At present there is no a-priori analysis to guide the selection */
/*     of NCV relative to NEV.  The only formal requrement is that NCV > NEV + 2. */
/*     However, it is recommended that NCV .ge. 2*NEV+1.  If many problems of */
/*     the same type are to be solved, one should experiment with increasing */
/*     NCV while keeping NEV fixed for a given test problem.  This will */
/*     usually decrease the required number of OP*x operations but it */
/*     also increases the work and storage required to maintain the orthogonal */
/*     basis vectors.  The optimal "cross-over" with respect to CPU time */
/*     is problem dependent and must be determined empirically. */
/*     See Chapter 8 of Reference 2 for further information. */

/*  5. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the */
/*     NP = IPARAM(8) real and imaginary parts of the shifts in locations */
/*         real part                  imaginary part */
/*         -----------------------    -------------- */
/*     1   WORKL(IPNTR(14))           WORKL(IPNTR(14)+NP) */
/*     2   WORKL(IPNTR(14)+1)         WORKL(IPNTR(14)+NP+1) */
/*                        .                          . */
/*                        .                          . */
/*                        .                          . */
/*     NP  WORKL(IPNTR(14)+NP-1)      WORKL(IPNTR(14)+2*NP-1). */

/*     Only complex conjugate pairs of shifts may be applied and the pairs */
/*     must be placed in consecutive locations. The real part of the */
/*     eigenvalues of the current upper Hessenberg matrix are located in */
/*     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1) and the imaginary part */
/*     in WORKL(IPNTR(7)) through WORKL(IPNTR(7)+NCV-1). They are ordered */
/*     according to the order defined by WHICH. The complex conjugate */
/*     pairs are kept together and the associated Ritz estimates are located in */
/*     WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1). */

/* ----------------------------------------------------------------------- */

/* \Data Distribution Note: */

/*  Fortran-D syntax: */
/*  ================ */
/*  Double precision resid(n), v(ldv,ncv), workd(3*n), workl(lworkl) */
/*  decompose  d1(n), d2(n,ncv) */
/*  align      resid(i) with d1(i) */
/*  align      v(i,j)   with d2(i,j) */
/*  align      workd(i) with d1(i)     range (1:n) */
/*  align      workd(i) with d1(i-n)   range (n+1:2*n) */
/*  align      workd(i) with d1(i-2*n) range (2*n+1:3*n) */
/*  distribute d1(block), d2(block,:) */
/*  replicated workl(lworkl) */

/*  Cray MPP syntax: */
/*  =============== */
/*  Double precision  resid(n), v(ldv,ncv), workd(n,3), workl(lworkl) */
/*  shared     resid(block), v(block,:), workd(block,:) */
/*  replicated workl(lworkl) */

/*  CM2/CM5 syntax: */
/*  ============== */

/* ----------------------------------------------------------------------- */

/*     include   'ex-nonsym.doc' */

/* ----------------------------------------------------------------------- */

/* \BeginLib */

/* \Local variables: */
/*     xxxxxx  real */

/* \References: */
/*  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in */
/*     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992), */
/*     pp 357-385. */
/*  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly */
/*     Restarted Arnoldi Iteration", Rice University Technical Report */
/*     TR95-13, Department of Computational and Applied Mathematics. */
/*  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for */
/*     Real Matrices", Linear Algebra and its Applications, vol 88/89, */
/*     pp 575-595, (1987). */

/* \Routines called: */
/*     dnaup2  ARPACK routine that implements the Implicitly Restarted */
/*             Arnoldi Iteration. */
/*     ivout   ARPACK utility routine that prints integers. */
/*     second  ARPACK utility routine for timing. */
/*     dvout   ARPACK utility routine that prints vectors. */
/*     dlamch  LAPACK routine that determines machine constants. */

/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \Revision history: */
/*     12/16/93: Version '1.1' */

/* \SCCS Information: @(#) */
/* FILE: naupd.F   SID: 2.10   DATE OF SID: 08/23/02   RELEASE: 2 */

/* \Remarks */

/* \EndLib */

/* ----------------------------------------------------------------------- */

int
dnaupd_ (integer * ido, char *bmat, integer * n, char *which, integer * nev, doublereal * tol, doublereal * resid, integer * ncv,
			doublereal * v, integer * ldv, integer * iparam, integer * ipntr,
			doublereal * workd, doublereal * workl, integer * lworkl, integer * info, ftnlen bmat_len, ftnlen which_len)
{
	/* Format strings */
	static char fmt_1000[] = "(//,5x,\002==================================="
		"==========\002,/5x,\002= Nonsymmetric implicit Arnoldi update co"
		"de =\002,/5x,\002= Version Number: \002,\002 2.4\002,21x,\002 "
		"=\002,/5x,\002= Version Date:   \002,\002 07/31/96\002,16x,\002 ="
		"\002,/5x,\002=============================================\002,/"
		"5x,\002= Summary of timing statistics              =\002,/5x," "\002=============================================\002,//)";
	static char fmt_1100[] = "(5x,\002Total number update iterations        "
		"     = \002,i5,/5x,\002Total number of OP*x operations          "
		"  = \002,i5,/5x,\002Total number of B*x operations             = "
		"\002,i5,/5x,\002Total number of reorthogonalization steps  = "
		"\002,i5,/5x,\002Total number of iterative refinement steps = "
		"\002,i5,/5x,\002Total number of restart steps              = "
		"\002,i5,/5x,\002Total time in user OP*x operation          = "
		"\002,f12.6,/5x,\002Total time in user B*x operation           ="
		" \002,f12.6,/5x,\002Total time in Arnoldi update routine       = "
		"\002,f12.6,/5x,\002Total time in naup2 routine                ="
		" \002,f12.6,/5x,\002Total time in basic Arnoldi iteration loop = "
		"\002,f12.6,/5x,\002Total time in reorthogonalization phase    ="
		" \002,f12.6,/5x,\002Total time in (re)start vector generation  = "
		"\002,f12.6,/5x,\002Total time in Hessenberg eig. subproblem   ="
		" \002,f12.6,/5x,\002Total time in getting the shifts           = "
		"\002,f12.6,/5x,\002Total time in applying the shifts          ="
		" \002,f12.6,/5x,\002Total time in convergence testing          = " "\002,f12.6,/5x,\002Total time in computing final Ritz vectors =" " \002,f12.6/)";

	/* System generated locals */
	integer v_dim1, v_offset, i__1, i__2;

	/* Builtin functions */
	integer s_cmp (char *, char *, ftnlen, ftnlen), s_wsfe (cilist *), e_wsfe (void), do_fio (integer *, char *, ftnlen);

	/* Local variables */
	static integer j;
	static real t0, t1;
	static integer nb, ih, iq, np, iw, ldh, ldq, nev0, mode, ierr, iupd, next, ritzi;
	static integer ritzr;
	static integer bounds, ishift, msglvl, mxiter;

	/* Fortran I/O blocks */
	static cilist io___1001 = { 0, 6, 0, fmt_1000, 0 };
	static cilist io___1002 = { 0, 6, 0, fmt_1100, 0 };



	/*     %----------------------------------------------------% */
	/*     | Include files for debugging and timing information | */
	/*     %----------------------------------------------------% */


	/* \SCCS Information: @(#) */
	/* FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 */

	/*     %---------------------------------% */
	/*     | See debug.doc for documentation | */
	/*     %---------------------------------% */

	/*     %------------------% */
	/*     | Scalar Arguments | */
	/*     %------------------% */

	/*     %--------------------------------% */
	/*     | See stat.doc for documentation | */
	/*     %--------------------------------% */

	/* \SCCS Information: @(#) */
	/* FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 */



	/*     %-----------------% */
	/*     | Array Arguments | */
	/*     %-----------------% */


	/*     %------------% */
	/*     | Parameters | */
	/*     %------------% */


	/*     %---------------% */
	/*     | Local Scalars | */
	/*     %---------------% */


	/*     %----------------------% */
	/*     | External Subroutines | */
	/*     %----------------------% */


	/*     %--------------------% */
	/*     | External Functions | */
	/*     %--------------------% */


	/*     %-----------------------% */
	/*     | Executable Statements | */
	/*     %-----------------------% */

	/* Parameter adjustments */
	--workd;
	--resid;
	v_dim1 = *ldv;
	v_offset = 1 + v_dim1;
	v -= v_offset;
	--iparam;
	--ipntr;
	--workl;

	/* Function Body */
	if (*ido == 0)
	{

		/*        %-------------------------------% */
		/*        | Initialize timing statistics  | */
		/*        | & message level for debugging | */
		/*        %-------------------------------% */

		dstatn_ ();
		second_ (&t0);
		msglvl = debug_1.mnaupd;

		/*        %----------------% */
		/*        | Error checking | */
		/*        %----------------% */

		ierr = 0;
		ishift = iparam[1];
		/*         levec  = iparam(2) */
		mxiter = iparam[3];
		/*         nb     = iparam(4) */
		nb = 1;

		/*        %--------------------------------------------% */
		/*        | Revision 2 performs only implicit restart. | */
		/*        %--------------------------------------------% */

		iupd = 1;
		mode = iparam[7];

		if (*n <= 0)
		{
			ierr = -1;
		}
		else if (*nev <= 0)
		{
			ierr = -2;
		}
		else if (*ncv <= *nev + 1 || *ncv > *n)
		{
			ierr = -3;
		}
		else if (mxiter <= 0)
		{
			ierr = 4;
		}
		else if (s_cmp (which, "LM", (ftnlen) 2, (ftnlen) 2) != 0 && s_cmp (which, "SM", (ftnlen) 2, (ftnlen) 2) != 0 && s_cmp (which, "LR",
																																										(ftnlen) 2, (ftnlen) 2) != 0
					&& s_cmp (which, "SR", (ftnlen) 2, (ftnlen) 2) != 0 && s_cmp (which, "LI", (ftnlen) 2, (ftnlen) 2) != 0
					&& s_cmp (which, "SI", (ftnlen) 2, (ftnlen) 2) != 0)
		{
			ierr = -5;
		}
		else if (*(unsigned char *) bmat != 'I' && *(unsigned char *) bmat != 'G')
		{
			ierr = -6;
		}
		else	/* if(complicated condition) */
		{
			/* Computing 2nd power */
			i__1 = *ncv;
			if (*lworkl < i__1 * i__1 * 3 + *ncv * 6)
			{
				ierr = -7;
			}
			else if (mode < 1 || mode > 4)
			{
				ierr = -10;
			}
			else if (mode == 1 && *(unsigned char *) bmat == 'G')
			{
				ierr = -11;
			}
			else if (ishift < 0 || ishift > 1)
			{
				ierr = -12;
			}
		}

		/*        %------------% */
		/*        | Error Exit | */
		/*        %------------% */

		if (ierr != 0)
		{
			*info = ierr;
			*ido = 99;
			goto L9000;
		}

		/*        %------------------------% */
		/*        | Set default parameters | */
		/*        %------------------------% */

		if (nb <= 0)
		{
			nb = 1;
		}
		if (*tol <= 0.)
		{
			*tol = dlamch_ ("EpsMach", (ftnlen) 7);
		}

		/*        %----------------------------------------------% */
		/*        | NP is the number of additional steps to      | */
		/*        | extend the length NEV Lanczos factorization. | */
		/*        | NEV0 is the local variable designating the   | */
		/*        | size of the invariant subspace desired.      | */
		/*        %----------------------------------------------% */

		np = *ncv - *nev;
		nev0 = *nev;

		/*        %-----------------------------% */
		/*        | Zero out internal workspace | */
		/*        %-----------------------------% */

		/* Computing 2nd power */
		i__2 = *ncv;
		i__1 = i__2 * i__2 * 3 + *ncv * 6;
		for (j = 1; j <= i__1; ++j)
		{
			workl[j] = 0.;
			/* L10: */
		}

		/*        %-------------------------------------------------------------% */
		/*        | Pointer into WORKL for address of H, RITZ, BOUNDS, Q        | */
		/*        | etc... and the remaining workspace.                         | */
		/*        | Also update pointer to be used on output.                   | */
		/*        | Memory is laid out as follows:                              | */
		/*        | workl(1:ncv*ncv) := generated Hessenberg matrix             | */
		/*        | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary        | */
		/*        |                                   parts of ritz values      | */
		/*        | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds        | */
		/*        | workl(ncv*ncv+3*ncv+1:2*ncv*ncv+3*ncv) := rotation matrix Q | */
		/*        | workl(2*ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) := workspace       | */
		/*        | The final workspace is needed by subroutine dneigh called   | */
		/*        | by dnaup2. Subroutine dneigh calls LAPACK routines for      | */
		/*        | calculating eigenvalues and the last row of the eigenvector | */
		/*        | matrix.                                                     | */
		/*        %-------------------------------------------------------------% */

		ldh = *ncv;
		ldq = *ncv;
		ih = 1;
		ritzr = ih + ldh * *ncv;
		ritzi = ritzr + *ncv;
		bounds = ritzi + *ncv;
		iq = bounds + *ncv;
		iw = iq + ldq * *ncv;
		/* Computing 2nd power */
		i__1 = *ncv;
		next = iw + i__1 * i__1 + *ncv * 3;

		ipntr[4] = next;
		ipntr[5] = ih;
		ipntr[6] = ritzr;
		ipntr[7] = ritzi;
		ipntr[8] = bounds;
		ipntr[14] = iw;

	}

	/*     %-------------------------------------------------------% */
	/*     | Carry out the Implicitly restarted Arnoldi Iteration. | */
	/*     %-------------------------------------------------------% */

	dnaup2_ (ido, bmat, n, which, &nev0, &np, tol, &resid[1], &mode, &iupd, &ishift, &mxiter, &v[v_offset], ldv, &workl[ih], &ldh, &workl[ritzr], &workl[ritzi],
				&workl[bounds], &workl[iq], &ldq, &workl[iw], &ipntr[1], &workd[1], info, (ftnlen) 1, (ftnlen) 2);

	/*     %--------------------------------------------------% */
	/*     | ido .ne. 99 implies use of reverse communication | */
	/*     | to compute operations involving OP or shifts.    | */
	/*     %--------------------------------------------------% */

	if (*ido == 3)
	{
		iparam[8] = np;
	}
	if (*ido != 99)
	{
		goto L9000;
	}

	iparam[3] = mxiter;
	iparam[5] = np;
	iparam[9] = timing_1.nopx;
	iparam[10] = timing_1.nbx;
	iparam[11] = timing_1.nrorth;

	/*     %------------------------------------% */
	/*     | Exit if there was an informational | */
	/*     | error within dnaup2.               | */
	/*     %------------------------------------% */

	if (*info < 0)
	{
		goto L9000;
	}
	if (*info == 2)
	{
		*info = 3;
	}

	if (msglvl > 0)
	{
		ivout_ (&debug_1.logfil, &c__1, &mxiter, &debug_1.ndigit, "_naupd: Nu" "mber of update iterations taken", (ftnlen) 41);
		ivout_ (&debug_1.logfil, &c__1, &np, &debug_1.ndigit, "_naupd: Number" " of wanted \"converged\" Ritz values", (ftnlen) 48);
		dvout_ (&debug_1.logfil, &np, &workl[ritzr], &debug_1.ndigit, "_naupd" ": Real part of the final Ritz values", (ftnlen) 42);
		dvout_ (&debug_1.logfil, &np, &workl[ritzi], &debug_1.ndigit, "_naupd" ": Imaginary part of the final Ritz values", (ftnlen) 47);
		dvout_ (&debug_1.logfil, &np, &workl[bounds], &debug_1.ndigit, "_naup" "d: Associated Ritz estimates", (ftnlen) 33);
	}

	second_ (&t1);
	timing_1.tnaupd = t1 - t0;

	if (msglvl > 0)
	{

		/*        %--------------------------------------------------------% */
		/*        | Version Number & Version Date are defined in version.h | */
		/*        %--------------------------------------------------------% */

		s_wsfe (&io___1001);
		e_wsfe ();
		s_wsfe (&io___1002);
		do_fio (&c__1, (char *) &mxiter, (ftnlen) sizeof (integer));
		do_fio (&c__1, (char *) &timing_1.nopx, (ftnlen) sizeof (integer));
		do_fio (&c__1, (char *) &timing_1.nbx, (ftnlen) sizeof (integer));
		do_fio (&c__1, (char *) &timing_1.nrorth, (ftnlen) sizeof (integer));
		do_fio (&c__1, (char *) &timing_1.nitref, (ftnlen) sizeof (integer));
		do_fio (&c__1, (char *) &timing_1.nrstrt, (ftnlen) sizeof (integer));
		do_fio (&c__1, (char *) &timing_1.tmvopx, (ftnlen) sizeof (real));
		do_fio (&c__1, (char *) &timing_1.tmvbx, (ftnlen) sizeof (real));
		do_fio (&c__1, (char *) &timing_1.tnaupd, (ftnlen) sizeof (real));
		do_fio (&c__1, (char *) &timing_1.tnaup2, (ftnlen) sizeof (real));
		do_fio (&c__1, (char *) &timing_1.tnaitr, (ftnlen) sizeof (real));
		do_fio (&c__1, (char *) &timing_1.titref, (ftnlen) sizeof (real));
		do_fio (&c__1, (char *) &timing_1.tgetv0, (ftnlen) sizeof (real));
		do_fio (&c__1, (char *) &timing_1.tneigh, (ftnlen) sizeof (real));
		do_fio (&c__1, (char *) &timing_1.tngets, (ftnlen) sizeof (real));
		do_fio (&c__1, (char *) &timing_1.tnapps, (ftnlen) sizeof (real));
		do_fio (&c__1, (char *) &timing_1.tnconv, (ftnlen) sizeof (real));
		do_fio (&c__1, (char *) &timing_1.trvec, (ftnlen) sizeof (real));
		e_wsfe ();
	}

 L9000:

	return 0;

	/*     %---------------% */
	/*     | End of dnaupd | */
	/*     %---------------% */

}	/* dnaupd_ */

/* ----------------------------------------------------------------------- */
/* \BeginDoc */

/* \Name: dnconv */

/* \Description: */
/*  Convergence testing for the nonsymmetric Arnoldi eigenvalue routine. */

/* \Usage: */
/*  call dnconv */
/*     ( N, RITZR, RITZI, BOUNDS, TOL, NCONV ) */

/* \Arguments */
/*  N       Integer.  (INPUT) */
/*          Number of Ritz values to check for convergence. */

/*  RITZR,  Double precision arrays of length N.  (INPUT) */
/*  RITZI   Real and imaginary parts of the Ritz values to be checked */
/*          for convergence. */
/*  BOUNDS  Double precision array of length N.  (INPUT) */
/*          Ritz estimates for the Ritz values in RITZR and RITZI. */

/*  TOL     Double precision scalar.  (INPUT) */
/*          Desired backward error for a Ritz value to be considered */
/*          "converged". */

/*  NCONV   Integer scalar.  (OUTPUT) */
/*          Number of "converged" Ritz values. */

/* \EndDoc */

/* ----------------------------------------------------------------------- */

/* \BeginLib */

/* \Local variables: */
/*     xxxxxx  real */

/* \Routines called: */
/*     second  ARPACK utility routine for timing. */
/*     dlamch  LAPACK routine that determines machine constants. */
/*     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully. */

/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \Revision history: */
/*     xx/xx/92: Version ' 2.1' */

/* \SCCS Information: @(#) */
/* FILE: nconv.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2 */

/* \Remarks */
/*     1. xxxx */

/* \EndLib */

/* ----------------------------------------------------------------------- */

static int
dnconv_ (integer * n, doublereal * ritzr, doublereal * ritzi, doublereal * bounds, doublereal * tol, integer * nconv)
{
	/* System generated locals */
	integer i__1;
	doublereal d__1, d__2;

	/* Builtin functions */
	double pow_dd (doublereal *, doublereal *);

	/* Local variables */
	static integer i__;
	static real t0, t1;
	static doublereal eps23, temp;


	/*     %----------------------------------------------------% */
	/*     | Include files for debugging and timing information | */
	/*     %----------------------------------------------------% */


	/* \SCCS Information: @(#) */
	/* FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 */

	/*     %---------------------------------% */
	/*     | See debug.doc for documentation | */
	/*     %---------------------------------% */

	/*     %------------------% */
	/*     | Scalar Arguments | */
	/*     %------------------% */

	/*     %--------------------------------% */
	/*     | See stat.doc for documentation | */
	/*     %--------------------------------% */

	/* \SCCS Information: @(#) */
	/* FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 */



	/*     %-----------------% */
	/*     | Array Arguments | */
	/*     %-----------------% */

	/*     %---------------% */
	/*     | Local Scalars | */
	/*     %---------------% */


	/*     %--------------------% */
	/*     | External Functions | */
	/*     %--------------------% */

	/*     %-----------------------% */
	/*     | Executable Statements | */
	/*     %-----------------------% */

	/*     %-------------------------------------------------------------% */
	/*     | Convergence test: unlike in the symmetric code, I am not    | */
	/*     | using things like refined error bounds and gap condition    | */
	/*     | because I don't know the exact equivalent concept.          | */
	/*     |                                                             | */
	/*     | Instead the i-th Ritz value is considered "converged" when: | */
	/*     |                                                             | */
	/*     |     bounds(i) .le. ( TOL * | ritz | )                       | */
	/*     |                                                             | */
	/*     | for some appropriate choice of norm.                        | */
	/*     %-------------------------------------------------------------% */

	/* Parameter adjustments */
	--bounds;
	--ritzi;
	--ritzr;

	/* Function Body */
	second_ (&t0);

	/*     %---------------------------------% */
	/*     | Get machine dependent constant. | */
	/*     %---------------------------------% */

	eps23 = dlamch_ ("Epsilon-Machine", (ftnlen) 15);
	eps23 = pow_dd (&eps23, &c_b2616);

	*nconv = 0;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		/* Computing MAX */
		d__1 = eps23, d__2 = dlapy2_ (&ritzr[i__], &ritzi[i__]);
		temp = max (d__1, d__2);
		if (bounds[i__] <= *tol * temp)
		{
			++(*nconv);
		}
		/* L20: */
	}

	second_ (&t1);
	timing_1.tnconv += t1 - t0;

	return 0;

	/*     %---------------% */
	/*     | End of dnconv | */
	/*     %---------------% */

}	/* dnconv_ */

/* ----------------------------------------------------------------------- */
/* \BeginDoc */

/* \Name: dneigh */

/* \Description: */
/*  Compute the eigenvalues of the current upper Hessenberg matrix */
/*  and the corresponding Ritz estimates given the current residual norm. */

/* \Usage: */
/*  call dneigh */
/*     ( RNORM, N, H, LDH, RITZR, RITZI, BOUNDS, Q, LDQ, WORKL, IERR ) */

/* \Arguments */
/*  RNORM   Double precision scalar.  (INPUT) */
/*          Residual norm corresponding to the current upper Hessenberg */
/*          matrix H. */

/*  N       Integer.  (INPUT) */
/*          Size of the matrix H. */

/*  H       Double precision N by N array.  (INPUT) */
/*          H contains the current upper Hessenberg matrix. */

/*  LDH     Integer.  (INPUT) */
/*          Leading dimension of H exactly as declared in the calling */
/*          program. */

/*  RITZR,  Double precision arrays of length N.  (OUTPUT) */
/*  RITZI   On output, RITZR(1:N) (resp. RITZI(1:N)) contains the real */
/*          (respectively imaginary) parts of the eigenvalues of H. */

/*  BOUNDS  Double precision array of length N.  (OUTPUT) */
/*          On output, BOUNDS contains the Ritz estimates associated with */
/*          the eigenvalues RITZR and RITZI.  This is equal to RNORM */
/*          times the last components of the eigenvectors corresponding */
/*          to the eigenvalues in RITZR and RITZI. */

/*  Q       Double precision N by N array.  (WORKSPACE) */
/*          Workspace needed to store the eigenvectors of H. */

/*  LDQ     Integer.  (INPUT) */
/*          Leading dimension of Q exactly as declared in the calling */
/*          program. */

/*  WORKL   Double precision work array of length N**2 + 3*N.  (WORKSPACE) */
/*          Private (replicated) array on each PE or array allocated on */
/*          the front end.  This is needed to keep the full Schur form */
/*          of H and also in the calculation of the eigenvectors of H. */

/*  IERR    Integer.  (OUTPUT) */
/*          Error exit flag from dlaqrb or dtrevc. */

/* \EndDoc */

/* ----------------------------------------------------------------------- */

/* \BeginLib */

/* \Local variables: */
/*     xxxxxx  real */

/* \Routines called: */
/*     dlaqrb  ARPACK routine to compute the real Schur form of an */
/*             upper Hessenberg matrix and last row of the Schur vectors. */
/*     second  ARPACK utility routine for timing. */
/*     dmout   ARPACK utility routine that prints matrices */
/*     dvout   ARPACK utility routine that prints vectors. */
/*     dlacpy  LAPACK matrix copy routine. */
/*     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully. */
/*     dtrevc  LAPACK routine to compute the eigenvectors of a matrix */
/*             in upper quasi-triangular form */
/*     dgemv   Level 2 BLAS routine for matrix vector multiplication. */
/*     dcopy   Level 1 BLAS that copies one vector to another . */
/*     dnrm2   Level 1 BLAS that computes the norm of a vector. */
/*     dscal   Level 1 BLAS that scales a vector. */


/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \Revision history: */
/*     xx/xx/92: Version ' 2.1' */

/* \SCCS Information: @(#) */
/* FILE: neigh.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2 */

/* \Remarks */
/*     None */

/* \EndLib */

/* ----------------------------------------------------------------------- */

static int
dneigh_ (doublereal * rnorm, integer * n, doublereal * h__,
			integer * ldh, doublereal * ritzr, doublereal * ritzi, doublereal * bounds, doublereal * q, integer * ldq, doublereal * workl, integer * ierr)
{
	/* System generated locals */
	integer h_dim1, h_offset, q_dim1, q_offset, i__1;
	doublereal d__1, d__2;

	/* Local variables */
	static integer i__;
	static real t0, t1;
	static doublereal vl[1], temp;
	static integer iconj;
	static logical select[1];
	static integer msglvl;

	/*     %----------------------------------------------------% */
	/*     | Include files for debugging and timing information | */
	/*     %----------------------------------------------------% */


	/* \SCCS Information: @(#) */
	/* FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 */

	/*     %---------------------------------% */
	/*     | See debug.doc for documentation | */
	/*     %---------------------------------% */

	/*     %------------------% */
	/*     | Scalar Arguments | */
	/*     %------------------% */

	/*     %--------------------------------% */
	/*     | See stat.doc for documentation | */
	/*     %--------------------------------% */

	/* \SCCS Information: @(#) */
	/* FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 */



	/*     %-----------------% */
	/*     | Array Arguments | */
	/*     %-----------------% */


	/*     %------------% */
	/*     | Parameters | */
	/*     %------------% */


	/*     %------------------------% */
	/*     | Local Scalars & Arrays | */
	/*     %------------------------% */


	/*     %----------------------% */
	/*     | External Subroutines | */
	/*     %----------------------% */


	/*     %--------------------% */
	/*     | External Functions | */
	/*     %--------------------% */


	/*     %---------------------% */
	/*     | Intrinsic Functions | */
	/*     %---------------------% */


	/*     %-----------------------% */
	/*     | Executable Statements | */
	/*     %-----------------------% */


	/*     %-------------------------------% */
	/*     | Initialize timing statistics  | */
	/*     | & message level for debugging | */
	/*     %-------------------------------% */

	/* Parameter adjustments */
	--workl;
	--bounds;
	--ritzi;
	--ritzr;
	h_dim1 = *ldh;
	h_offset = 1 + h_dim1;
	h__ -= h_offset;
	q_dim1 = *ldq;
	q_offset = 1 + q_dim1;
	q -= q_offset;

	/* Function Body */
	second_ (&t0);
	msglvl = debug_1.mneigh;

	if (msglvl > 2)
	{
		dmout_ (&debug_1.logfil, n, n, &h__[h_offset], ldh, &debug_1.ndigit, "_neigh: Entering upper Hessenberg matrix H ", (ftnlen) 43);
	}

	/*     %-----------------------------------------------------------% */
	/*     | 1. Compute the eigenvalues, the last components of the    | */
	/*     |    corresponding Schur vectors and the full Schur form T  | */
	/*     |    of the current upper Hessenberg matrix H.              | */
	/*     | dlaqrb returns the full Schur form of H in WORKL(1:N**2)  | */
	/*     | and the last components of the Schur vectors in BOUNDS.   | */
	/*     %-----------------------------------------------------------% */

	dlacpy_ ("All", n, n, &h__[h_offset], ldh, &workl[1], n, (ftnlen) 3);
	dlaqrb_ (&c_true, n, &c__1, n, &workl[1], n, &ritzr[1], &ritzi[1], &bounds[1], ierr);
	if (*ierr != 0)
	{
		goto L9000;
	}

	if (msglvl > 1)
	{
		dvout_ (&debug_1.logfil, n, &bounds[1], &debug_1.ndigit, "_neigh: las" "t row of the Schur matrix for H", (ftnlen) 42);
	}

	/*     %-----------------------------------------------------------% */
	/*     | 2. Compute the eigenvectors of the full Schur form T and  | */
	/*     |    apply the last components of the Schur vectors to get  | */
	/*     |    the last components of the corresponding eigenvectors. | */
	/*     | Remember that if the i-th and (i+1)-st eigenvalues are    | */
	/*     | complex conjugate pairs, then the real & imaginary part   | */
	/*     | of the eigenvector components are split across adjacent   | */
	/*     | columns of Q.                                             | */
	/*     %-----------------------------------------------------------% */

	dtrevc_ ("R", "A", select, n, &workl[1], n, vl, n, &q[q_offset], ldq, n, n, &workl[*n * *n + 1], ierr, (ftnlen) 1, (ftnlen) 1);

	if (*ierr != 0)
	{
		goto L9000;
	}

	/*     %------------------------------------------------% */
	/*     | Scale the returning eigenvectors so that their | */
	/*     | euclidean norms are all one. LAPACK subroutine | */
	/*     | dtrevc returns each eigenvector normalized so  | */
	/*     | that the element of largest magnitude has      | */
	/*     | magnitude 1; here the magnitude of a complex   | */
	/*     | number (x,y) is taken to be |x| + |y|.         | */
	/*     %------------------------------------------------% */

	iconj = 0;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		if ((d__1 = ritzi[i__], abs (d__1)) <= 0.)
		{

			/*           %----------------------% */
			/*           | Real eigenvalue case | */
			/*           %----------------------% */

			temp = dnrm2_ (n, &q[i__ * q_dim1 + 1], &c__1);
			d__1 = 1. / temp;
			dscal_ (n, &d__1, &q[i__ * q_dim1 + 1], &c__1);
		}
		else
		{

			/*           %-------------------------------------------% */
			/*           | Complex conjugate pair case. Note that    | */
			/*           | since the real and imaginary part of      | */
			/*           | the eigenvector are stored in consecutive | */
			/*           | columns, we further normalize by the      | */
			/*           | square root of two.                       | */
			/*           %-------------------------------------------% */

			if (iconj == 0)
			{
				d__1 = dnrm2_ (n, &q[i__ * q_dim1 + 1], &c__1);
				d__2 = dnrm2_ (n, &q[(i__ + 1) * q_dim1 + 1], &c__1);
				temp = dlapy2_ (&d__1, &d__2);
				d__1 = 1. / temp;
				dscal_ (n, &d__1, &q[i__ * q_dim1 + 1], &c__1);
				d__1 = 1. / temp;
				dscal_ (n, &d__1, &q[(i__ + 1) * q_dim1 + 1], &c__1);
				iconj = 1;
			}
			else
			{
				iconj = 0;
			}
		}
		/* L10: */
	}

	dgemv_ ("T", n, n, &c_b348, &q[q_offset], ldq, &bounds[1], &c__1, &c_b507, &workl[1], &c__1, (ftnlen) 1);

	if (msglvl > 1)
	{
		dvout_ (&debug_1.logfil, n, &workl[1], &debug_1.ndigit, "_neigh: Last" " row of the eigenvector matrix for H", (ftnlen) 48);
	}

	/*     %----------------------------% */
	/*     | Compute the Ritz estimates | */
	/*     %----------------------------% */

	iconj = 0;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		if ((d__1 = ritzi[i__], abs (d__1)) <= 0.)
		{

			/*           %----------------------% */
			/*           | Real eigenvalue case | */
			/*           %----------------------% */

			bounds[i__] = *rnorm * (d__1 = workl[i__], abs (d__1));
		}
		else
		{

			/*           %-------------------------------------------% */
			/*           | Complex conjugate pair case. Note that    | */
			/*           | since the real and imaginary part of      | */
			/*           | the eigenvector are stored in consecutive | */
			/*           | columns, we need to take the magnitude    | */
			/*           | of the last components of the two vectors | */
			/*           %-------------------------------------------% */

			if (iconj == 0)
			{
				bounds[i__] = *rnorm * dlapy2_ (&workl[i__], &workl[i__ + 1]);
				bounds[i__ + 1] = bounds[i__];
				iconj = 1;
			}
			else
			{
				iconj = 0;
			}
		}
		/* L20: */
	}

	if (msglvl > 2)
	{
		dvout_ (&debug_1.logfil, n, &ritzr[1], &debug_1.ndigit, "_neigh: Real" " part of the eigenvalues of H", (ftnlen) 41);
		dvout_ (&debug_1.logfil, n, &ritzi[1], &debug_1.ndigit, "_neigh: Imag" "inary part of the eigenvalues of H", (ftnlen) 46);
		dvout_ (&debug_1.logfil, n, &bounds[1], &debug_1.ndigit, "_neigh: Rit" "z estimates for the eigenvalues of H", (ftnlen) 47);
	}

	second_ (&t1);
	timing_1.tneigh += t1 - t0;

 L9000:
	return 0;

	/*     %---------------% */
	/*     | End of dneigh | */
	/*     %---------------% */

}	/* dneigh_ */

/* \BeginDoc */

/* \Name: dneupd */

/* \Description: */

/*  This subroutine returns the converged approximations to eigenvalues */
/*  of A*z = lambda*B*z and (optionally): */

/*      (1) The corresponding approximate eigenvectors; */

/*      (2) An orthonormal basis for the associated approximate */
/*          invariant subspace; */

/*      (3) Both. */

/*  There is negligible additional cost to obtain eigenvectors.  An orthonormal */
/*  basis is always computed.  There is an additional storage cost of n*nev */
/*  if both are requested (in this case a separate array Z must be supplied). */

/*  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z */
/*  are derived from approximate eigenvalues and eigenvectors of */
/*  of the linear operator OP prescribed by the MODE selection in the */
/*  call to DNAUPD .  DNAUPD  must be called before this routine is called. */
/*  These approximate eigenvalues and vectors are commonly called Ritz */
/*  values and Ritz vectors respectively.  They are referred to as such */
/*  in the comments that follow.  The computed orthonormal basis for the */
/*  invariant subspace corresponding to these Ritz values is referred to as a */
/*  Schur basis. */

/*  See documentation in the header of the subroutine DNAUPD  for */
/*  definition of OP as well as other terms and the relation of computed */
/*  Ritz values and Ritz vectors of OP with respect to the given problem */
/*  A*z = lambda*B*z.  For a brief description, see definitions of */
/*  IPARAM(7), MODE and WHICH in the documentation of DNAUPD . */

/* \Usage: */
/*  call dneupd */
/*     ( RVEC, HOWMNY, SELECT, DR, DI, Z, LDZ, SIGMAR, SIGMAI, WORKEV, BMAT, */
/*       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, */
/*       LWORKL, INFO ) */

/* \Arguments: */
/*  RVEC    LOGICAL  (INPUT) */
/*          Specifies whether a basis for the invariant subspace corresponding */
/*          to the converged Ritz value approximations for the eigenproblem */
/*          A*z = lambda*B*z is computed. */

/*             RVEC = .FALSE.     Compute Ritz values only. */

/*             RVEC = .TRUE.      Compute the Ritz vectors or Schur vectors. */
/*                                See Remarks below. */

/*  HOWMNY  Character*1  (INPUT) */
/*          Specifies the form of the basis for the invariant subspace */
/*          corresponding to the converged Ritz values that is to be computed. */

/*          = 'A': Compute NEV Ritz vectors; */
/*          = 'P': Compute NEV Schur vectors; */
/*          = 'S': compute some of the Ritz vectors, specified */
/*                 by the logical array SELECT. */

/*  SELECT  Logical array of dimension NCV.  (INPUT) */
/*          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be */
/*          computed. To select the Ritz vector corresponding to a */
/*          Ritz value (DR(j), DI(j)), SELECT(j) must be set to .TRUE.. */
/*          If HOWMNY = 'A' or 'P', SELECT is used as internal workspace. */

/*  DR      Double precision  array of dimension NEV+1.  (OUTPUT) */
/*          If IPARAM(7) = 1,2 or 3 and SIGMAI=0.0  then on exit: DR contains */
/*          the real part of the Ritz  approximations to the eigenvalues of */
/*          A*z = lambda*B*z. */
/*          If IPARAM(7) = 3, 4 and SIGMAI is not equal to zero, then on exit: */
/*          DR contains the real part of the Ritz values of OP computed by */
/*          DNAUPD . A further computation must be performed by the user */
/*          to transform the Ritz values computed for OP by DNAUPD  to those */
/*          of the original system A*z = lambda*B*z. See remark 3 below. */

/*  DI      Double precision  array of dimension NEV+1.  (OUTPUT) */
/*          On exit, DI contains the imaginary part of the Ritz value */
/*          approximations to the eigenvalues of A*z = lambda*B*z associated */
/*          with DR. */

/*          NOTE: When Ritz values are complex, they will come in complex */
/*                conjugate pairs.  If eigenvectors are requested, the */
/*                corresponding Ritz vectors will also come in conjugate */
/*                pairs and the real and imaginary parts of these are */
/*                represented in two consecutive columns of the array Z */
/*                (see below). */

/*  Z       Double precision  N by NEV+1 array if RVEC = .TRUE. and HOWMNY = 'A'. (OUTPUT) */
/*          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of */
/*          Z represent approximate eigenvectors (Ritz vectors) corresponding */
/*          to the NCONV=IPARAM(5) Ritz values for eigensystem */
/*          A*z = lambda*B*z. */

/*          The complex Ritz vector associated with the Ritz value */
/*          with positive imaginary part is stored in two consecutive */
/*          columns.  The first column holds the real part of the Ritz */
/*          vector and the second column holds the imaginary part.  The */
/*          Ritz vector associated with the Ritz value with negative */
/*          imaginary part is simply the complex conjugate of the Ritz vector */
/*          associated with the positive imaginary part. */

/*          If  RVEC = .FALSE. or HOWMNY = 'P', then Z is not referenced. */

/*          NOTE: If if RVEC = .TRUE. and a Schur basis is not required, */
/*          the array Z may be set equal to first NEV+1 columns of the Arnoldi */
/*          basis array V computed by DNAUPD .  In this case the Arnoldi basis */
/*          will be destroyed and overwritten with the eigenvector basis. */

/*  LDZ     Integer.  (INPUT) */
/*          The leading dimension of the array Z.  If Ritz vectors are */
/*          desired, then  LDZ >= max( 1, N ).  In any case,  LDZ >= 1. */

/*  SIGMAR  Double precision   (INPUT) */
/*          If IPARAM(7) = 3 or 4, represents the real part of the shift. */
/*          Not referenced if IPARAM(7) = 1 or 2. */

/*  SIGMAI  Double precision   (INPUT) */
/*          If IPARAM(7) = 3 or 4, represents the imaginary part of the shift. */
/*          Not referenced if IPARAM(7) = 1 or 2. See remark 3 below. */

/*  WORKEV  Double precision  work array of dimension 3*NCV.  (WORKSPACE) */

/*  **** The remaining arguments MUST be the same as for the   **** */
/*  **** call to DNAUPD  that was just completed.               **** */

/*  NOTE: The remaining arguments */

/*           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, */
/*           WORKD, WORKL, LWORKL, INFO */

/*         must be passed directly to DNEUPD  following the last call */
/*         to DNAUPD .  These arguments MUST NOT BE MODIFIED between */
/*         the the last call to DNAUPD  and the call to DNEUPD . */

/*  Three of these parameters (V, WORKL, INFO) are also output parameters: */

/*  V       Double precision  N by NCV array.  (INPUT/OUTPUT) */

/*          Upon INPUT: the NCV columns of V contain the Arnoldi basis */
/*                      vectors for OP as constructed by DNAUPD  . */

/*          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns */
/*                       contain approximate Schur vectors that span the */
/*                       desired invariant subspace.  See Remark 2 below. */

/*          NOTE: If the array Z has been set equal to first NEV+1 columns */
/*          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the */
/*          Arnoldi basis held by V has been overwritten by the desired */
/*          Ritz vectors.  If a separate array Z has been passed then */
/*          the first NCONV=IPARAM(5) columns of V will contain approximate */
/*          Schur vectors that span the desired invariant subspace. */

/*  WORKL   Double precision  work array of length LWORKL.  (OUTPUT/WORKSPACE) */
/*          WORKL(1:ncv*ncv+3*ncv) contains information obtained in */
/*          dnaupd .  They are not changed by dneupd . */
/*          WORKL(ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) holds the */
/*          real and imaginary part of the untransformed Ritz values, */
/*          the upper quasi-triangular matrix for H, and the */
/*          associated matrix representation of the invariant subspace for H. */

/*          Note: IPNTR(9:13) contains the pointer into WORKL for addresses */
/*          of the above information computed by dneupd . */
/*          ------------------------------------------------------------- */
/*          IPNTR(9):  pointer to the real part of the NCV RITZ values of the */
/*                     original system. */
/*          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of */
/*                     the original system. */
/*          IPNTR(11): pointer to the NCV corresponding error bounds. */
/*          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular */
/*                     Schur matrix for H. */
/*          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors */
/*                     of the upper Hessenberg matrix H. Only referenced by */
/*                     dneupd  if RVEC = .TRUE. See Remark 2 below. */
/*          ------------------------------------------------------------- */

/*  INFO    Integer.  (OUTPUT) */
/*          Error flag on output. */

/*          =  0: Normal exit. */

/*          =  1: The Schur form computed by LAPACK routine dlahqr */
/*                could not be reordered by LAPACK routine dtrsen . */
/*                Re-enter subroutine dneupd  with IPARAM(5)=NCV and */
/*                increase the size of the arrays DR and DI to have */
/*                dimension at least dimension NCV and allocate at least NCV */
/*                columns for Z. NOTE: Not necessary if Z and V share */
/*                the same space. Please notify the authors if this error */
/*                occurs. */

/*          = -1: N must be positive. */
/*          = -2: NEV must be positive. */
/*          = -3: NCV-NEV >= 2 and less than or equal to N. */
/*          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI' */
/*          = -6: BMAT must be one of 'I' or 'G'. */
/*          = -7: Length of private work WORKL array is not sufficient. */
/*          = -8: Error return from calculation of a real Schur form. */
/*                Informational error from LAPACK routine dlahqr . */
/*          = -9: Error return from calculation of eigenvectors. */
/*                Informational error from LAPACK routine dtrevc . */
/*          = -10: IPARAM(7) must be 1,2,3,4. */
/*          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible. */
/*          = -12: HOWMNY = 'S' not yet implemented */
/*          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true. */
/*          = -14: DNAUPD  did not find any eigenvalues to sufficient */
/*                 accuracy. */
/*          = -15: DNEUPD got a different count of the number of converged */
/*                 Ritz values than DNAUPD got.  This indicates the user */
/*                 probably made an error in passing data from DNAUPD to */
/*                 DNEUPD or that the data was modified before entering */
/*                 DNEUPD */

/* \BeginLib */

/* \References: */
/*  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in */
/*     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992), */
/*     pp 357-385. */
/*  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly */
/*     Restarted Arnoldi Iteration", Rice University Technical Report */
/*     TR95-13, Department of Computational and Applied Mathematics. */
/*  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for */
/*     Real Matrices", Linear Algebra and its Applications, vol 88/89, */
/*     pp 575-595, (1987). */

/* \Routines called: */
/*     ivout   ARPACK utility routine that prints integers. */
/*     dmout    ARPACK utility routine that prints matrices */
/*     dvout    ARPACK utility routine that prints vectors. */
/*     dgeqr2   LAPACK routine that computes the QR factorization of */
/*             a matrix. */
/*     dlacpy   LAPACK matrix copy routine. */
/*     dlahqr   LAPACK routine to compute the real Schur form of an */
/*             upper Hessenberg matrix. */
/*     dlamch   LAPACK routine that determines machine constants. */
/*     dlapy2   LAPACK routine to compute sqrt(x**2+y**2) carefully. */
/*     dlaset   LAPACK matrix initialization routine. */
/*     dorm2r   LAPACK routine that applies an orthogonal matrix in */
/*             factored form. */
/*     dtrevc   LAPACK routine to compute the eigenvectors of a matrix */
/*             in upper quasi-triangular form. */
/*     dtrsen   LAPACK routine that re-orders the Schur form. */
/*     dtrmm    Level 3 BLAS matrix times an upper triangular matrix. */
/*     dger     Level 2 BLAS rank one update to a matrix. */
/*     dcopy    Level 1 BLAS that copies one vector to another . */
/*     ddot     Level 1 BLAS that computes the scalar product of two vectors. */
/*     dnrm2    Level 1 BLAS that computes the norm of a vector. */
/*     dscal    Level 1 BLAS that scales a vector. */

/* \Remarks */

/*  1. Currently only HOWMNY = 'A' and 'P' are implemented. */

/*     Let trans(X) denote the transpose of X. */

/*  2. Schur vectors are an orthogonal representation for the basis of */
/*     Ritz vectors. Thus, their numerical properties are often superior. */
/*     If RVEC = .TRUE. then the relationship */
/*             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and */
/*     trans(V(:,1:IPARAM(5))) * V(:,1:IPARAM(5)) = I are approximately */
/*     satisfied. Here T is the leading submatrix of order IPARAM(5) of the */
/*     real upper quasi-triangular matrix stored workl(ipntr(12)). That is, */
/*     T is block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; */
/*     each 2-by-2 diagonal block has its diagonal elements equal and its */
/*     off-diagonal elements of opposite sign.  Corresponding to each 2-by-2 */
/*     diagonal block is a complex conjugate pair of Ritz values. The real */
/*     Ritz values are stored on the diagonal of T. */

/*  3. If IPARAM(7) = 3 or 4 and SIGMAI is not equal zero, then the user must */
/*     form the IPARAM(5) Rayleigh quotients in order to transform the Ritz */
/*     values computed by DNAUPD  for OP to those of A*z = lambda*B*z. */
/*     Set RVEC = .true. and HOWMNY = 'A', and */
/*     compute */
/*           trans(Z(:,I)) * A * Z(:,I) if DI(I) = 0. */
/*     If DI(I) is not equal to zero and DI(I+1) = - D(I), */
/*     then the desired real and imaginary parts of the Ritz value are */
/*           trans(Z(:,I)) * A * Z(:,I) +  trans(Z(:,I+1)) * A * Z(:,I+1), */
/*           trans(Z(:,I)) * A * Z(:,I+1) -  trans(Z(:,I+1)) * A * Z(:,I), */
/*     respectively. */
/*     Another possibility is to set RVEC = .true. and HOWMNY = 'P' and */
/*     compute trans(V(:,1:IPARAM(5))) * A * V(:,1:IPARAM(5)) and then an upper */
/*     quasi-triangular matrix of order IPARAM(5) is computed. See remark */
/*     2 above. */

/* \Authors */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Chao Yang                    Houston, Texas */
/*     Dept. of Computational & */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \SCCS Information: @(#) */
/* FILE: neupd.F   SID: 2.7   DATE OF SID: 09/20/00   RELEASE: 2 */

/* \EndLib */

/* ----------------------------------------------------------------------- */
int
dneupd_ (logical * rvec, char *howmny, logical * select,
			doublereal * dr, doublereal * di, doublereal * z__, integer * ldz,
			doublereal * sigmar, doublereal * sigmai, doublereal * workev, char *bmat, integer * n, char *which, integer * nev, doublereal * tol,
			doublereal * resid, integer * ncv, doublereal * v, integer * ldv, integer
			* iparam, integer * ipntr, doublereal * workd, doublereal * workl,
			integer * lworkl, integer * info, ftnlen howmny_len, ftnlen bmat_len, ftnlen which_len)
{
	/* System generated locals */
	integer v_dim1, v_offset, z_dim1, z_offset, i__1;
	doublereal d__1, d__2;

	/* Builtin functions */
	double pow_dd (doublereal *, doublereal *);
	integer s_cmp (char *, char *, ftnlen, ftnlen);
	/* Subroutine */ int s_copy (char *, char *, ftnlen, ftnlen);

	/* Local variables */
	static integer j, k, ih, jj, np;
	static doublereal vl[1] /* was [1][1] */ ;
	static integer ibd, ldh, ldq, iri;
	static doublereal sep;
	static integer irr, wri, wrr;
	static integer mode;
	static doublereal eps23;
	static integer ierr;
	static doublereal temp;
	static integer iwev;
	static char type__[6];
	static doublereal temp1;
	static integer ihbds, iconj;
	static doublereal conds;
	static logical reord;
	static integer nconv;
	static integer iwork[1];
	static doublereal rnorm;
	static integer ritzi;
	static integer ritzr;
	static integer iheigi, iheigr, bounds, invsub, iuptri, msglvl, outncv, ishift, numcnv;


	/*     %----------------------------------------------------% */
	/*     | Include files for debugging and timing information | */
	/*     %----------------------------------------------------% */


	/* \SCCS Information: @(#) */
	/* FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 */

	/*     %---------------------------------% */
	/*     | See debug.doc for documentation | */
	/*     %---------------------------------% */

	/*     %------------------% */
	/*     | Scalar Arguments | */
	/*     %------------------% */

	/*     %--------------------------------% */
	/*     | See stat.doc for documentation | */
	/*     %--------------------------------% */

	/* \SCCS Information: @(#) */
	/* FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 */



	/*     %-----------------% */
	/*     | Array Arguments | */
	/*     %-----------------% */


	/*     %------------% */
	/*     | Parameters | */
	/*     %------------% */


	/*     %---------------% */
	/*     | Local Scalars | */
	/*     %---------------% */


	/*     %----------------------% */
	/*     | External Subroutines | */
	/*     %----------------------% */


	/*     %--------------------% */
	/*     | External Functions | */
	/*     %--------------------% */


	/*     %---------------------% */
	/*     | Intrinsic Functions | */
	/*     %---------------------% */


	/*     %-----------------------% */
	/*     | Executable Statements | */
	/*     %-----------------------% */

	/*     %------------------------% */
	/*     | Set default parameters | */
	/*     %------------------------% */

	/* Parameter adjustments */
	z_dim1 = *ldz;
	z_offset = 1 + z_dim1;
	z__ -= z_offset;
	--workd;
	--resid;
	--di;
	--dr;
	--workev;
	--select;
	v_dim1 = *ldv;
	v_offset = 1 + v_dim1;
	v -= v_offset;
	--iparam;
	--ipntr;
	--workl;

	/* Function Body */
	msglvl = debug_1.mneupd;
	mode = iparam[7];
	nconv = iparam[5];
	*info = 0;

	/*     %---------------------------------% */
	/*     | Get machine dependent constant. | */
	/*     %---------------------------------% */

	eps23 = dlamch_ ("Epsilon-Machine", (ftnlen) 15);
	eps23 = pow_dd (&eps23, &c_b2616);

	/*     %--------------% */
	/*     | Quick return | */
	/*     %--------------% */

	ierr = 0;

	if (nconv <= 0)
	{
		ierr = -14;
	}
	else if (*n <= 0)
	{
		ierr = -1;
	}
	else if (*nev <= 0)
	{
		ierr = -2;
	}
	else if (*ncv <= *nev + 1 || *ncv > *n)
	{
		ierr = -3;
	}
	else if (s_cmp (which, "LM", (ftnlen) 2, (ftnlen) 2) != 0 && s_cmp (which,
																							  "SM", (ftnlen) 2, (ftnlen) 2) != 0 && s_cmp (which, "LR", (ftnlen) 2,
																																						  (ftnlen) 2) != 0
				&& s_cmp (which, "SR", (ftnlen) 2, (ftnlen) 2) != 0 && s_cmp (which, "LI", (ftnlen) 2, (ftnlen) 2) != 0
				&& s_cmp (which, "SI", (ftnlen) 2, (ftnlen) 2) != 0)
	{
		ierr = -5;
	}
	else if (*(unsigned char *) bmat != 'I' && *(unsigned char *) bmat != 'G')
	{
		ierr = -6;
	}
	else	/* if(complicated condition) */
	{
		/* Computing 2nd power */
		i__1 = *ncv;
		if (*lworkl < i__1 * i__1 * 3 + *ncv * 6)
		{
			ierr = -7;
		}
		else if (*(unsigned char *) howmny != 'A' && *(unsigned char *) howmny != 'P' && *(unsigned char *) howmny != 'S' && *rvec)
		{
			ierr = -13;
		}
		else if (*(unsigned char *) howmny == 'S')
		{
			ierr = -12;
		}
	}

	if (mode == 1 || mode == 2)
	{
		s_copy (type__, "REGULR", (ftnlen) 6, (ftnlen) 6);
	}
	else if (mode == 3 && *sigmai == 0.)
	{
		s_copy (type__, "SHIFTI", (ftnlen) 6, (ftnlen) 6);
	}
	else if (mode == 3)
	{
		s_copy (type__, "REALPT", (ftnlen) 6, (ftnlen) 6);
	}
	else if (mode == 4)
	{
		s_copy (type__, "IMAGPT", (ftnlen) 6, (ftnlen) 6);
	}
	else
	{
		ierr = -10;
	}
	if (mode == 1 && *(unsigned char *) bmat == 'G')
	{
		ierr = -11;
	}

	/*     %------------% */
	/*     | Error Exit | */
	/*     %------------% */

	if (ierr != 0)
	{
		*info = ierr;
		goto L9000;
	}

	/*     %--------------------------------------------------------% */
	/*     | Pointer into WORKL for address of H, RITZ, BOUNDS, Q   | */
	/*     | etc... and the remaining workspace.                    | */
	/*     | Also update pointer to be used on output.              | */
	/*     | Memory is laid out as follows:                         | */
	/*     | workl(1:ncv*ncv) := generated Hessenberg matrix        | */
	/*     | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary   | */
	/*     |                                   parts of ritz values | */
	/*     | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds   | */
	/*     %--------------------------------------------------------% */

	/*     %-----------------------------------------------------------% */
	/*     | The following is used and set by DNEUPD .                  | */
	/*     | workl(ncv*ncv+3*ncv+1:ncv*ncv+4*ncv) := The untransformed | */
	/*     |                             real part of the Ritz values. | */
	/*     | workl(ncv*ncv+4*ncv+1:ncv*ncv+5*ncv) := The untransformed | */
	/*     |                        imaginary part of the Ritz values. | */
	/*     | workl(ncv*ncv+5*ncv+1:ncv*ncv+6*ncv) := The untransformed | */
	/*     |                           error bounds of the Ritz values | */
	/*     | workl(ncv*ncv+6*ncv+1:2*ncv*ncv+6*ncv) := Holds the upper | */
	/*     |                             quasi-triangular matrix for H | */
	/*     | workl(2*ncv*ncv+6*ncv+1: 3*ncv*ncv+6*ncv) := Holds the    | */
	/*     |       associated matrix representation of the invariant   | */
	/*     |       subspace for H.                                     | */
	/*     | GRAND total of NCV * ( 3 * NCV + 6 ) locations.           | */
	/*     %-----------------------------------------------------------% */

	ih = ipntr[5];
	ritzr = ipntr[6];
	ritzi = ipntr[7];
	bounds = ipntr[8];
	ldh = *ncv;
	ldq = *ncv;
	iheigr = bounds + ldh;
	iheigi = iheigr + ldh;
	ihbds = iheigi + ldh;
	iuptri = ihbds + ldh;
	invsub = iuptri + ldh * *ncv;
	ipntr[9] = iheigr;
	ipntr[10] = iheigi;
	ipntr[11] = ihbds;
	ipntr[12] = iuptri;
	ipntr[13] = invsub;
	wrr = 1;
	wri = *ncv + 1;
	iwev = wri + *ncv;

	/*     %-----------------------------------------% */
	/*     | irr points to the REAL part of the Ritz | */
	/*     |     values computed by _neigh before    | */
	/*     |     exiting _naup2.                     | */
	/*     | iri points to the IMAGINARY part of the | */
	/*     |     Ritz values computed by _neigh      | */
	/*     |     before exiting _naup2.              | */
	/*     | ibd points to the Ritz estimates        | */
	/*     |     computed by _neigh before exiting   | */
	/*     |     _naup2.                             | */
	/*     %-----------------------------------------% */

	irr = ipntr[14] + *ncv * *ncv;
	iri = irr + *ncv;
	ibd = iri + *ncv;

	/*     %------------------------------------% */
	/*     | RNORM is B-norm of the RESID(1:N). | */
	/*     %------------------------------------% */

	rnorm = workl[ih + 2];
	workl[ih + 2] = 0.;

	if (msglvl > 2)
	{
		dvout_ (&debug_1.logfil, ncv, &workl[irr], &debug_1.ndigit, "_neupd: " "Real part of Ritz values passed in from _NAUPD.", (ftnlen) 55);
		dvout_ (&debug_1.logfil, ncv, &workl[iri], &debug_1.ndigit, "_neupd: " "Imag part of Ritz values passed in from _NAUPD.", (ftnlen) 55);
		dvout_ (&debug_1.logfil, ncv, &workl[ibd], &debug_1.ndigit, "_neupd: " "Ritz estimates passed in from _NAUPD.", (ftnlen) 45);
	}

	if (*rvec)
	{

		reord = FALSE_;

		/*        %---------------------------------------------------% */
		/*        | Use the temporary bounds array to store indices   | */
		/*        | These will be used to mark the select array later | */
		/*        %---------------------------------------------------% */

		i__1 = *ncv;
		for (j = 1; j <= i__1; ++j)
		{
			workl[bounds + j - 1] = (doublereal) j;
			select[j] = FALSE_;
			/* L10: */
		}

		/*        %-------------------------------------% */
		/*        | Select the wanted Ritz values.      | */
		/*        | Sort the Ritz values so that the    | */
		/*        | wanted ones appear at the tailing   | */
		/*        | NEV positions of workl(irr) and     | */
		/*        | workl(iri).  Move the corresponding | */
		/*        | error estimates in workl(bound)     | */
		/*        | accordingly.                        | */
		/*        %-------------------------------------% */

		np = *ncv - *nev;
		ishift = 0;
		dngets_ (&ishift, which, nev, &np, &workl[irr], &workl[iri], &workl[bounds], &workl[1], &workl[np + 1], (ftnlen) 2);

		if (msglvl > 2)
		{
			dvout_ (&debug_1.logfil, ncv, &workl[irr], &debug_1.ndigit, "_neu" "pd: Real part of Ritz values after calling _NGETS.", (ftnlen) 54);
			dvout_ (&debug_1.logfil, ncv, &workl[iri], &debug_1.ndigit, "_neu" "pd: Imag part of Ritz values after calling _NGETS.", (ftnlen) 54);
			dvout_ (&debug_1.logfil, ncv, &workl[bounds], &debug_1.ndigit, "_neupd: Ritz value indices after calling _NGETS.", (ftnlen) 48);
		}

		/*        %-----------------------------------------------------% */
		/*        | Record indices of the converged wanted Ritz values  | */
		/*        | Mark the select array for possible reordering       | */
		/*        %-----------------------------------------------------% */

		numcnv = 0;
		i__1 = *ncv;
		for (j = 1; j <= i__1; ++j)
		{
			/* Computing MAX */
			d__1 = eps23, d__2 = dlapy2_ (&workl[irr + *ncv - j], &workl[iri + *ncv - j]);
			temp1 = max (d__1, d__2);
			jj = (integer) workl[bounds + *ncv - j];
			if (numcnv < nconv && workl[ibd + jj - 1] <= *tol * temp1)
			{
				select[jj] = TRUE_;
				++numcnv;
				if (jj > *nev)
				{
					reord = TRUE_;
				}
			}
			/* L11: */
		}

		/*        %-----------------------------------------------------------% */
		/*        | Check the count (numcnv) of converged Ritz values with    | */
		/*        | the number (nconv) reported by dnaupd.  If these two      | */
		/*        | are different then there has probably been an error       | */
		/*        | caused by incorrect passing of the dnaupd data.           | */
		/*        %-----------------------------------------------------------% */

		if (msglvl > 2)
		{
			ivout_ (&debug_1.logfil, &c__1, &numcnv, &debug_1.ndigit, "_neupd" ": Number of specified eigenvalues", (ftnlen) 39);
			ivout_ (&debug_1.logfil, &c__1, &nconv, &debug_1.ndigit, "_neupd:" " Number of \"converged\" eigenvalues", (ftnlen) 41);
		}

		if (numcnv != nconv)
		{
			*info = -15;
			goto L9000;
		}

		/*        %-----------------------------------------------------------% */
		/*        | Call LAPACK routine dlahqr  to compute the real Schur form | */
		/*        | of the upper Hessenberg matrix returned by DNAUPD .        | */
		/*        | Make a copy of the upper Hessenberg matrix.               | */
		/*        | Initialize the Schur vector matrix Q to the identity.     | */
		/*        %-----------------------------------------------------------% */

		i__1 = ldh * *ncv;
		dcopy_ (&i__1, &workl[ih], &c__1, &workl[iuptri], &c__1);
		dlaset_ ("All", ncv, ncv, &c_b507, &c_b348, &workl[invsub], &ldq, (ftnlen) 3);
		dlahqr_ (&c_true, &c_true, ncv, &c__1, ncv, &workl[iuptri], &ldh, &workl[iheigr], &workl[iheigi], &c__1, ncv, &workl[invsub], &ldq, &ierr);
		dcopy_ (ncv, &workl[invsub + *ncv - 1], &ldq, &workl[ihbds], &c__1);

		if (ierr != 0)
		{
			*info = -8;
			goto L9000;
		}

		if (msglvl > 1)
		{
			dvout_ (&debug_1.logfil, ncv, &workl[iheigr], &debug_1.ndigit, "_neupd: Real part of the eigenvalues of H", (ftnlen) 41);
			dvout_ (&debug_1.logfil, ncv, &workl[iheigi], &debug_1.ndigit, "_neupd: Imaginary part of the Eigenvalues of H", (ftnlen) 46);
			dvout_ (&debug_1.logfil, ncv, &workl[ihbds], &debug_1.ndigit, "_neupd: Last row of the Schur vector matrix", (ftnlen) 43);
			if (msglvl > 3)
			{
				dmout_ (&debug_1.logfil, ncv, ncv, &workl[iuptri], &ldh, &debug_1.ndigit, "_neupd: The upper quasi-triangular " "matrix ", (ftnlen) 42);
			}
		}

		if (reord)
		{

			/*           %-----------------------------------------------------% */
			/*           | Reorder the computed upper quasi-triangular matrix. | */
			/*           %-----------------------------------------------------% */

			dtrsen_ ("None", "V", &select[1], ncv, &workl[iuptri], &ldh, &workl[invsub], &ldq, &workl[iheigr], &workl[iheigi], &nconv, &conds, &sep, &workl[ihbds],
						ncv, iwork, &c__1, &ierr, (ftnlen) 4, (ftnlen) 1);

			if (ierr == 1)
			{
				*info = 1;
				goto L9000;
			}

			if (msglvl > 2)
			{
				dvout_ (&debug_1.logfil, ncv, &workl[iheigr], &debug_1.ndigit, "_neupd: Real part of the eigenvalues of H--reordered", (ftnlen) 52);
				dvout_ (&debug_1.logfil, ncv, &workl[iheigi], &debug_1.ndigit, "_neupd: Imag part of the eigenvalues of H--reordered", (ftnlen) 52);
				if (msglvl > 3)
				{
					dmout_ (&debug_1.logfil, ncv, ncv, &workl[iuptri], &ldq, &debug_1.ndigit, "_neupd: Quasi-triangular matrix" " after re-ordering", (ftnlen) 49);
				}
			}

		}

		/*        %---------------------------------------% */
		/*        | Copy the last row of the Schur vector | */
		/*        | into workl(ihbds).  This will be used | */
		/*        | to compute the Ritz estimates of      | */
		/*        | converged Ritz values.                | */
		/*        %---------------------------------------% */

		dcopy_ (ncv, &workl[invsub + *ncv - 1], &ldq, &workl[ihbds], &c__1);

		/*        %----------------------------------------------------% */
		/*        | Place the computed eigenvalues of H into DR and DI | */
		/*        | if a spectral transformation was not used.         | */
		/*        %----------------------------------------------------% */

		if (s_cmp (type__, "REGULR", (ftnlen) 6, (ftnlen) 6) == 0)
		{
			dcopy_ (&nconv, &workl[iheigr], &c__1, &dr[1], &c__1);
			dcopy_ (&nconv, &workl[iheigi], &c__1, &di[1], &c__1);
		}

		/*        %----------------------------------------------------------% */
		/*        | Compute the QR factorization of the matrix representing  | */
		/*        | the wanted invariant subspace located in the first NCONV | */
		/*        | columns of workl(invsub,ldq).                            | */
		/*        %----------------------------------------------------------% */

		dgeqr2_ (ncv, &nconv, &workl[invsub], &ldq, &workev[1], &workev[*ncv + 1], &ierr);

		/*        %---------------------------------------------------------% */
		/*        | * Postmultiply V by Q using dorm2r .                     | */
		/*        | * Copy the first NCONV columns of VQ into Z.            | */
		/*        | * Postmultiply Z by R.                                  | */
		/*        | The N by NCONV matrix Z is now a matrix representation  | */
		/*        | of the approximate invariant subspace associated with   | */
		/*        | the Ritz values in workl(iheigr) and workl(iheigi)      | */
		/*        | The first NCONV columns of V are now approximate Schur  | */
		/*        | vectors associated with the real upper quasi-triangular | */
		/*        | matrix of order NCONV in workl(iuptri)                  | */
		/*        %---------------------------------------------------------% */

		dorm2r_ ("Right", "Notranspose", n, ncv, &nconv, &workl[invsub], &ldq, &workev[1], &v[v_offset], ldv, &workd[*n + 1], &ierr, (ftnlen) 5, (ftnlen) 11);
		dlacpy_ ("All", n, &nconv, &v[v_offset], ldv, &z__[z_offset], ldz, (ftnlen) 3);

		i__1 = nconv;
		for (j = 1; j <= i__1; ++j)
		{

			/*           %---------------------------------------------------% */
			/*           | Perform both a column and row scaling if the      | */
			/*           | diagonal element of workl(invsub,ldq) is negative | */
			/*           | I'm lazy and don't take advantage of the upper    | */
			/*           | quasi-triangular form of workl(iuptri,ldq)        | */
			/*           | Note that since Q is orthogonal, R is a diagonal  | */
			/*           | matrix consisting of plus or minus ones           | */
			/*           %---------------------------------------------------% */

			if (workl[invsub + (j - 1) * ldq + j - 1] < 0.)
			{
				dscal_ (&nconv, &c_b347, &workl[iuptri + j - 1], &ldq);
				dscal_ (&nconv, &c_b347, &workl[iuptri + (j - 1) * ldq], &c__1);
			}

			/* L20: */
		}

		if (*(unsigned char *) howmny == 'A')
		{

			/*           %--------------------------------------------% */
			/*           | Compute the NCONV wanted eigenvectors of T | */
			/*           | located in workl(iuptri,ldq).              | */
			/*           %--------------------------------------------% */

			i__1 = *ncv;
			for (j = 1; j <= i__1; ++j)
			{
				if (j <= nconv)
				{
					select[j] = TRUE_;
				}
				else
				{
					select[j] = FALSE_;
				}
				/* L30: */
			}

			dtrevc_ ("Right", "Select", &select[1], ncv, &workl[iuptri], &ldq,
						vl, &c__1, &workl[invsub], &ldq, ncv, &outncv, &workev[1], &ierr, (ftnlen) 5, (ftnlen) 6);

			if (ierr != 0)
			{
				*info = -9;
				goto L9000;
			}

			/*           %------------------------------------------------% */
			/*           | Scale the returning eigenvectors so that their | */
			/*           | Euclidean norms are all one. LAPACK subroutine | */
			/*           | dtrevc  returns each eigenvector normalized so  | */
			/*           | that the element of largest magnitude has      | */
			/*           | magnitude 1;                                   | */
			/*           %------------------------------------------------% */

			iconj = 0;
			i__1 = nconv;
			for (j = 1; j <= i__1; ++j)
			{

				if (workl[iheigi + j - 1] == 0.)
				{

					/*                 %----------------------% */
					/*                 | real eigenvalue case | */
					/*                 %----------------------% */

					temp = dnrm2_ (ncv, &workl[invsub + (j - 1) * ldq], &c__1);
					d__1 = 1. / temp;
					dscal_ (ncv, &d__1, &workl[invsub + (j - 1) * ldq], &c__1);

				}
				else
				{

					/*                 %-------------------------------------------% */
					/*                 | Complex conjugate pair case. Note that    | */
					/*                 | since the real and imaginary part of      | */
					/*                 | the eigenvector are stored in consecutive | */
					/*                 | columns, we further normalize by the      | */
					/*                 | square root of two.                       | */
					/*                 %-------------------------------------------% */

					if (iconj == 0)
					{
						d__1 = dnrm2_ (ncv, &workl[invsub + (j - 1) * ldq], &c__1);
						d__2 = dnrm2_ (ncv, &workl[invsub + j * ldq], &c__1);
						temp = dlapy2_ (&d__1, &d__2);
						d__1 = 1. / temp;
						dscal_ (ncv, &d__1, &workl[invsub + (j - 1) * ldq], &c__1);
						d__1 = 1. / temp;
						dscal_ (ncv, &d__1, &workl[invsub + j * ldq], &c__1);
						iconj = 1;
					}
					else
					{
						iconj = 0;
					}

				}

				/* L40: */
			}

			dgemv_ ("T", ncv, &nconv, &c_b348, &workl[invsub], &ldq, &workl[ihbds], &c__1, &c_b507, &workev[1], &c__1, (ftnlen) 1);

			iconj = 0;
			i__1 = nconv;
			for (j = 1; j <= i__1; ++j)
			{
				if (workl[iheigi + j - 1] != 0.)
				{

					/*                 %-------------------------------------------% */
					/*                 | Complex conjugate pair case. Note that    | */
					/*                 | since the real and imaginary part of      | */
					/*                 | the eigenvector are stored in consecutive | */
					/*                 %-------------------------------------------% */

					if (iconj == 0)
					{
						workev[j] = dlapy2_ (&workev[j], &workev[j + 1]);
						workev[j + 1] = workev[j];
						iconj = 1;
					}
					else
					{
						iconj = 0;
					}
				}
				/* L45: */
			}

			if (msglvl > 2)
			{
				dcopy_ (ncv, &workl[invsub + *ncv - 1], &ldq, &workl[ihbds], &c__1);
				dvout_ (&debug_1.logfil, ncv, &workl[ihbds], &debug_1.ndigit, "_neupd: Last row of the eigenvector matrix for T", (ftnlen) 48);
				if (msglvl > 3)
				{
					dmout_ (&debug_1.logfil, ncv, ncv, &workl[invsub], &ldq, &debug_1.ndigit, "_neupd: The eigenvector matrix " "for T", (ftnlen) 36);
				}
			}

			/*           %---------------------------------------% */
			/*           | Copy Ritz estimates into workl(ihbds) | */
			/*           %---------------------------------------% */

			dcopy_ (&nconv, &workev[1], &c__1, &workl[ihbds], &c__1);

			/*           %---------------------------------------------------------% */
			/*           | Compute the QR factorization of the eigenvector matrix  | */
			/*           | associated with leading portion of T in the first NCONV | */
			/*           | columns of workl(invsub,ldq).                           | */
			/*           %---------------------------------------------------------% */

			dgeqr2_ (ncv, &nconv, &workl[invsub], &ldq, &workev[1], &workev[*ncv + 1], &ierr);

			/*           %----------------------------------------------% */
			/*           | * Postmultiply Z by Q.                       | */
			/*           | * Postmultiply Z by R.                       | */
			/*           | The N by NCONV matrix Z is now contains the  | */
			/*           | Ritz vectors associated with the Ritz values | */
			/*           | in workl(iheigr) and workl(iheigi).          | */
			/*           %----------------------------------------------% */

			dorm2r_ ("Right", "Notranspose", n, ncv, &nconv, &workl[invsub], &ldq, &workev[1], &z__[z_offset], ldz, &workd[*n + 1], &ierr, (ftnlen) 5,
						(ftnlen) 11);

			dtrmm_ ("Right", "Upper", "No transpose", "Non-unit", n, &nconv, &c_b348, &workl[invsub], &ldq, &z__[z_offset], ldz, (ftnlen) 5, (ftnlen) 5,
					  (ftnlen) 12, (ftnlen) 8);

		}

	}
	else
	{

		/*        %------------------------------------------------------% */
		/*        | An approximate invariant subspace is not needed.     | */
		/*        | Place the Ritz values computed DNAUPD  into DR and DI | */
		/*        %------------------------------------------------------% */

		dcopy_ (&nconv, &workl[ritzr], &c__1, &dr[1], &c__1);
		dcopy_ (&nconv, &workl[ritzi], &c__1, &di[1], &c__1);
		dcopy_ (&nconv, &workl[ritzr], &c__1, &workl[iheigr], &c__1);
		dcopy_ (&nconv, &workl[ritzi], &c__1, &workl[iheigi], &c__1);
		dcopy_ (&nconv, &workl[bounds], &c__1, &workl[ihbds], &c__1);
	}

	/*     %------------------------------------------------% */
	/*     | Transform the Ritz values and possibly vectors | */
	/*     | and corresponding error bounds of OP to those  | */
	/*     | of A*x = lambda*B*x.                           | */
	/*     %------------------------------------------------% */

	if (s_cmp (type__, "REGULR", (ftnlen) 6, (ftnlen) 6) == 0)
	{

		if (*rvec)
		{
			dscal_ (ncv, &rnorm, &workl[ihbds], &c__1);
		}

	}
	else
	{

		/*        %---------------------------------------% */
		/*        |   A spectral transformation was used. | */
		/*        | * Determine the Ritz estimates of the | */
		/*        |   Ritz values in the original system. | */
		/*        %---------------------------------------% */

		if (s_cmp (type__, "SHIFTI", (ftnlen) 6, (ftnlen) 6) == 0)
		{

			if (*rvec)
			{
				dscal_ (ncv, &rnorm, &workl[ihbds], &c__1);
			}

			i__1 = *ncv;
			for (k = 1; k <= i__1; ++k)
			{
				temp = dlapy2_ (&workl[iheigr + k - 1], &workl[iheigi + k - 1]);
				workl[ihbds + k - 1] = (d__1 = workl[ihbds + k - 1], abs (d__1)) / temp / temp;
				/* L50: */
			}

		}
		else if (s_cmp (type__, "REALPT", (ftnlen) 6, (ftnlen) 6) == 0)
		{

			i__1 = *ncv;
			for (k = 1; k <= i__1; ++k)
			{
				/* L60: */
			}

		}
		else if (s_cmp (type__, "IMAGPT", (ftnlen) 6, (ftnlen) 6) == 0)
		{

			i__1 = *ncv;
			for (k = 1; k <= i__1; ++k)
			{
				/* L70: */
			}

		}

		/*        %-----------------------------------------------------------% */
		/*        | *  Transform the Ritz values back to the original system. | */
		/*        |    For TYPE = 'SHIFTI' the transformation is              | */
		/*        |             lambda = 1/theta + sigma                      | */
		/*        |    For TYPE = 'REALPT' or 'IMAGPT' the user must from     | */
		/*        |    Rayleigh quotients or a projection. See remark 3 above.| */
		/*        | NOTES:                                                    | */
		/*        | *The Ritz vectors are not affected by the transformation. | */
		/*        %-----------------------------------------------------------% */

		if (s_cmp (type__, "SHIFTI", (ftnlen) 6, (ftnlen) 6) == 0)
		{

			i__1 = *ncv;
			for (k = 1; k <= i__1; ++k)
			{
				temp = dlapy2_ (&workl[iheigr + k - 1], &workl[iheigi + k - 1]);
				workl[iheigr + k - 1] = workl[iheigr + k - 1] / temp / temp + *sigmar;
				workl[iheigi + k - 1] = -workl[iheigi + k - 1] / temp / temp + *sigmai;
				/* L80: */
			}

			dcopy_ (&nconv, &workl[iheigr], &c__1, &dr[1], &c__1);
			dcopy_ (&nconv, &workl[iheigi], &c__1, &di[1], &c__1);

		}
		else if (s_cmp (type__, "REALPT", (ftnlen) 6, (ftnlen) 6) == 0 || s_cmp (type__, "IMAGPT", (ftnlen) 6, (ftnlen) 6) == 0)
		{

			dcopy_ (&nconv, &workl[iheigr], &c__1, &dr[1], &c__1);
			dcopy_ (&nconv, &workl[iheigi], &c__1, &di[1], &c__1);

		}

	}

	if (s_cmp (type__, "SHIFTI", (ftnlen) 6, (ftnlen) 6) == 0 && msglvl > 1)
	{
		dvout_ (&debug_1.logfil, &nconv, &dr[1], &debug_1.ndigit, "_neupd: Un" "transformed real part of the Ritz valuess.", (ftnlen) 52);
		dvout_ (&debug_1.logfil, &nconv, &di[1], &debug_1.ndigit, "_neupd: Un" "transformed imag part of the Ritz valuess.", (ftnlen) 52);
		dvout_ (&debug_1.logfil, &nconv, &workl[ihbds], &debug_1.ndigit, "_ne" "upd: Ritz estimates of untransformed Ritz values.", (ftnlen) 52);
	}
	else if (s_cmp (type__, "REGULR", (ftnlen) 6, (ftnlen) 6) == 0 && msglvl > 1)
	{
		dvout_ (&debug_1.logfil, &nconv, &dr[1], &debug_1.ndigit, "_neupd: Re" "al parts of converged Ritz values.", (ftnlen) 44);
		dvout_ (&debug_1.logfil, &nconv, &di[1], &debug_1.ndigit, "_neupd: Im" "ag parts of converged Ritz values.", (ftnlen) 44);
		dvout_ (&debug_1.logfil, &nconv, &workl[ihbds], &debug_1.ndigit, "_ne" "upd: Associated Ritz estimates.", (ftnlen) 34);
	}

	/*     %-------------------------------------------------% */
	/*     | Eigenvector Purification step. Formally perform | */
	/*     | one of inverse subspace iteration. Only used    | */
	/*     | for MODE = 2.                                   | */
	/*     %-------------------------------------------------% */

	if (*rvec && *(unsigned char *) howmny == 'A' && s_cmp (type__, "SHIFTI", (ftnlen) 6, (ftnlen) 6) == 0)
	{

		/*        %------------------------------------------------% */
		/*        | Purify the computed Ritz vectors by adding a   | */
		/*        | little bit of the residual vector:             | */
		/*        |                      T                         | */
		/*        |          resid(:)*( e    s ) / theta           | */
		/*        |                      NCV                       | */
		/*        | where H s = s theta. Remember that when theta  | */
		/*        | has nonzero imaginary part, the corresponding  | */
		/*        | Ritz vector is stored across two columns of Z. | */
		/*        %------------------------------------------------% */

		iconj = 0;
		i__1 = nconv;
		for (j = 1; j <= i__1; ++j)
		{
			if (workl[iheigi + j - 1] == 0.)
			{
				workev[j] = workl[invsub + (j - 1) * ldq + *ncv - 1] / workl[iheigr + j - 1];
			}
			else if (iconj == 0)
			{
				temp = dlapy2_ (&workl[iheigr + j - 1], &workl[iheigi + j - 1]);
				workev[j] = (workl[invsub + (j - 1) * ldq + *ncv - 1] * workl[iheigr + j - 1] + workl[invsub + j * ldq + *ncv - 1] *
								 workl[iheigi + j - 1]) / temp / temp;
				workev[j + 1] = (workl[invsub + j * ldq + *ncv - 1] * workl[iheigr + j - 1] - workl[invsub + (j - 1) * ldq + *ncv
																																- 1] * workl[iheigi + j - 1]) / temp / temp;
				iconj = 1;
			}
			else
			{
				iconj = 0;
			}
			/* L110: */
		}

		/*        %---------------------------------------% */
		/*        | Perform a rank one update to Z and    | */
		/*        | purify all the Ritz vectors together. | */
		/*        %---------------------------------------% */

		dger_ (n, &nconv, &c_b348, &resid[1], &c__1, &workev[1], &c__1, &z__[z_offset], ldz);

	}

 L9000:

	return 0;

	/*     %---------------% */
	/*     | End of DNEUPD  | */
	/*     %---------------% */

}	/* dneupd_ */

/* ----------------------------------------------------------------------- */
/* \BeginDoc */

/* \Name: dngets */

/* \Description: */
/*  Given the eigenvalues of the upper Hessenberg matrix H, */
/*  computes the NP shifts AMU that are zeros of the polynomial of */
/*  degree NP which filters out components of the unwanted eigenvectors */
/*  corresponding to the AMU's based on some given criteria. */

/*  NOTE: call this even in the case of user specified shifts in order */
/*  to sort the eigenvalues, and error bounds of H for later use. */

/* \Usage: */
/*  call dngets */
/*     ( ISHIFT, WHICH, KEV, NP, RITZR, RITZI, BOUNDS, SHIFTR, SHIFTI ) */

/* \Arguments */
/*  ISHIFT  Integer.  (INPUT) */
/*          Method for selecting the implicit shifts at each iteration. */
/*          ISHIFT = 0: user specified shifts */
/*          ISHIFT = 1: exact shift with respect to the matrix H. */

/*  WHICH   Character*2.  (INPUT) */
/*          Shift selection criteria. */
/*          'LM' -> want the KEV eigenvalues of largest magnitude. */
/*          'SM' -> want the KEV eigenvalues of smallest magnitude. */
/*          'LR' -> want the KEV eigenvalues of largest real part. */
/*          'SR' -> want the KEV eigenvalues of smallest real part. */
/*          'LI' -> want the KEV eigenvalues of largest imaginary part. */
/*          'SI' -> want the KEV eigenvalues of smallest imaginary part. */

/*  KEV      Integer.  (INPUT/OUTPUT) */
/*           INPUT: KEV+NP is the size of the matrix H. */
/*           OUTPUT: Possibly increases KEV by one to keep complex conjugate */
/*           pairs together. */

/*  NP       Integer.  (INPUT/OUTPUT) */
/*           Number of implicit shifts to be computed. */
/*           OUTPUT: Possibly decreases NP by one to keep complex conjugate */
/*           pairs together. */

/*  RITZR,  Double precision array of length KEV+NP.  (INPUT/OUTPUT) */
/*  RITZI   On INPUT, RITZR and RITZI contain the real and imaginary */
/*          parts of the eigenvalues of H. */
/*          On OUTPUT, RITZR and RITZI are sorted so that the unwanted */
/*          eigenvalues are in the first NP locations and the wanted */
/*          portion is in the last KEV locations.  When exact shifts are */
/*          selected, the unwanted part corresponds to the shifts to */
/*          be applied. Also, if ISHIFT .eq. 1, the unwanted eigenvalues */
/*          are further sorted so that the ones with largest Ritz values */
/*          are first. */

/*  BOUNDS  Double precision array of length KEV+NP.  (INPUT/OUTPUT) */
/*          Error bounds corresponding to the ordering in RITZ. */

/*  SHIFTR, SHIFTI  *** USE deprecated as of version 2.1. *** */


/* \EndDoc */

/* ----------------------------------------------------------------------- */

/* \BeginLib */

/* \Local variables: */
/*     xxxxxx  real */

/* \Routines called: */
/*     dsortc  ARPACK sorting routine. */
/*     dcopy   Level 1 BLAS that copies one vector to another . */

/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \Revision history: */
/*     xx/xx/92: Version ' 2.1' */

/* \SCCS Information: @(#) */
/* FILE: ngets.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2 */

/* \Remarks */
/*     1. xxxx */

/* \EndLib */

/* ----------------------------------------------------------------------- */

static int
dngets_ (integer * ishift, char *which, integer * kev,
			integer * np, doublereal * ritzr, doublereal * ritzi, doublereal * bounds, doublereal * shiftr, doublereal * shifti, ftnlen which_len)
{
	/* System generated locals */
	integer i__1;

	/* Builtin functions */
	integer s_cmp (char *, char *, ftnlen, ftnlen);

	/* Local variables */
	static real t0, t1;
	static integer msglvl;


	/*     %----------------------------------------------------% */
	/*     | Include files for debugging and timing information | */
	/*     %----------------------------------------------------% */


	/* \SCCS Information: @(#) */
	/* FILE: debug.h   SID: 2.3   DATE OF SID: 11/16/95   RELEASE: 2 */

	/*     %---------------------------------% */
	/*     | See debug.doc for documentation | */
	/*     %---------------------------------% */

	/*     %------------------% */
	/*     | Scalar Arguments | */
	/*     %------------------% */

	/*     %--------------------------------% */
	/*     | See stat.doc for documentation | */
	/*     %--------------------------------% */

	/* \SCCS Information: @(#) */
	/* FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 */



	/*     %-----------------% */
	/*     | Array Arguments | */
	/*     %-----------------% */


	/*     %------------% */
	/*     | Parameters | */
	/*     %------------% */


	/*     %---------------% */
	/*     | Local Scalars | */
	/*     %---------------% */


	/*     %----------------------% */
	/*     | External Subroutines | */
	/*     %----------------------% */


	/*     %----------------------% */
	/*     | Intrinsics Functions | */
	/*     %----------------------% */


	/*     %-----------------------% */
	/*     | Executable Statements | */
	/*     %-----------------------% */

	/*     %-------------------------------% */
	/*     | Initialize timing statistics  | */
	/*     | & message level for debugging | */
	/*     %-------------------------------% */

	/* Parameter adjustments */
	--bounds;
	--ritzi;
	--ritzr;
	--shiftr;
	--shifti;

	/* Function Body */
	second_ (&t0);
	msglvl = debug_1.mngets;

	/*     %----------------------------------------------------% */
	/*     | LM, SM, LR, SR, LI, SI case.                       | */
	/*     | Sort the eigenvalues of H into the desired order   | */
	/*     | and apply the resulting order to BOUNDS.           | */
	/*     | The eigenvalues are sorted so that the wanted part | */
	/*     | are always in the last KEV locations.              | */
	/*     | We first do a pre-processing sort in order to keep | */
	/*     | complex conjugate pairs together                   | */
	/*     %----------------------------------------------------% */

	if (s_cmp (which, "LM", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		i__1 = *kev + *np;
		dsortc_ ("LR", &c_true, &i__1, &ritzr[1], &ritzi[1], &bounds[1], (ftnlen) 2);
	}
	else if (s_cmp (which, "SM", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		i__1 = *kev + *np;
		dsortc_ ("SR", &c_true, &i__1, &ritzr[1], &ritzi[1], &bounds[1], (ftnlen) 2);
	}
	else if (s_cmp (which, "LR", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		i__1 = *kev + *np;
		dsortc_ ("LM", &c_true, &i__1, &ritzr[1], &ritzi[1], &bounds[1], (ftnlen) 2);
	}
	else if (s_cmp (which, "SR", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		i__1 = *kev + *np;
		dsortc_ ("SM", &c_true, &i__1, &ritzr[1], &ritzi[1], &bounds[1], (ftnlen) 2);
	}
	else if (s_cmp (which, "LI", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		i__1 = *kev + *np;
		dsortc_ ("LM", &c_true, &i__1, &ritzr[1], &ritzi[1], &bounds[1], (ftnlen) 2);
	}
	else if (s_cmp (which, "SI", (ftnlen) 2, (ftnlen) 2) == 0)
	{
		i__1 = *kev + *np;
		dsortc_ ("SM", &c_true, &i__1, &ritzr[1], &ritzi[1], &bounds[1], (ftnlen) 2);
	}

	i__1 = *kev + *np;
	dsortc_ (which, &c_true, &i__1, &ritzr[1], &ritzi[1], &bounds[1], (ftnlen) 2);

	/*     %-------------------------------------------------------% */
	/*     | Increase KEV by one if the ( ritzr(np),ritzi(np) )    | */
	/*     | = ( ritzr(np+1),-ritzi(np+1) ) and ritz(np) .ne. zero | */
	/*     | Accordingly decrease NP by one. In other words keep   | */
	/*     | complex conjugate pairs together.                     | */
	/*     %-------------------------------------------------------% */

	if (ritzr[*np + 1] - ritzr[*np] == 0. && ritzi[*np + 1] + ritzi[*np] == 0.)
	{
		--(*np);
		++(*kev);
	}

	if (*ishift == 1)
	{

		/*        %-------------------------------------------------------% */
		/*        | Sort the unwanted Ritz values used as shifts so that  | */
		/*        | the ones with largest Ritz estimates are first        | */
		/*        | This will tend to minimize the effects of the         | */
		/*        | forward instability of the iteration when they shifts | */
		/*        | are applied in subroutine dnapps.                     | */
		/*        | Be careful and use 'SR' since we want to sort BOUNDS! | */
		/*        %-------------------------------------------------------% */

		dsortc_ ("SR", &c_true, np, &bounds[1], &ritzr[1], &ritzi[1], (ftnlen) 2);
	}

	second_ (&t1);
	timing_1.tngets += t1 - t0;

	if (msglvl > 0)
	{
		ivout_ (&debug_1.logfil, &c__1, kev, &debug_1.ndigit, "_ngets: KEV is", (ftnlen) 14);
		ivout_ (&debug_1.logfil, &c__1, np, &debug_1.ndigit, "_ngets: NP is", (ftnlen) 13);
		i__1 = *kev + *np;
		dvout_ (&debug_1.logfil, &i__1, &ritzr[1], &debug_1.ndigit, "_ngets: " "Eigenvalues of current H matrix -- real part", (ftnlen) 52);
		i__1 = *kev + *np;
		dvout_ (&debug_1.logfil, &i__1, &ritzi[1], &debug_1.ndigit, "_ngets: " "Eigenvalues of current H matrix -- imag part", (ftnlen) 52);
		i__1 = *kev + *np;
		dvout_ (&debug_1.logfil, &i__1, &bounds[1], &debug_1.ndigit, "_ngets:" " Ritz estimates of the current KEV+NP Ritz values", (ftnlen) 56);
	}

	return 0;

	/*     %---------------% */
	/*     | End of dngets | */
	/*     %---------------% */

}	/* dngets_ */

/* ----------------------------------------------------------------------- */
/* \BeginDoc */

/* \Name: dsortc */

/* \Description: */
/*  Sorts the complex array in XREAL and XIMAG into the order */
/*  specified by WHICH and optionally applies the permutation to the */
/*  real array Y. It is assumed that if an element of XIMAG is */
/*  nonzero, then its negative is also an element. In other words, */
/*  both members of a complex conjugate pair are to be sorted and the */
/*  pairs are kept adjacent to each other. */

/* \Usage: */
/*  call dsortc */
/*     ( WHICH, APPLY, N, XREAL, XIMAG, Y ) */

/* \Arguments */
/*  WHICH   Character*2.  (Input) */
/*          'LM' -> sort XREAL,XIMAG into increasing order of magnitude. */
/*          'SM' -> sort XREAL,XIMAG into decreasing order of magnitude. */
/*          'LR' -> sort XREAL into increasing order of algebraic. */
/*          'SR' -> sort XREAL into decreasing order of algebraic. */
/*          'LI' -> sort XIMAG into increasing order of magnitude. */
/*          'SI' -> sort XIMAG into decreasing order of magnitude. */
/*          NOTE: If an element of XIMAG is non-zero, then its negative */
/*                is also an element. */

/*  APPLY   Logical.  (Input) */
/*          APPLY = .TRUE.  -> apply the sorted order to array Y. */
/*          APPLY = .FALSE. -> do not apply the sorted order to array Y. */

/*  N       Integer.  (INPUT) */
/*          Size of the arrays. */

/*  XREAL,  Double precision array of length N.  (INPUT/OUTPUT) */
/*  XIMAG   Real and imaginary part of the array to be sorted. */

/*  Y       Double precision array of length N.  (INPUT/OUTPUT) */

/* \EndDoc */

/* ----------------------------------------------------------------------- */

/* \BeginLib */

/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \Revision history: */
/*     xx/xx/92: Version ' 2.1' */
/*               Adapted from the sort routine in LANSO. */

/* \SCCS Information: @(#) */
/* FILE: sortc.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2 */

/* \EndLib */

/* ----------------------------------------------------------------------- */

static int
dsortc_ (char *which, logical * apply, integer * n, doublereal * xreal, doublereal * ximag, doublereal * y, ftnlen which_len)
{
	/* System generated locals */
	integer i__1;
	doublereal d__1, d__2;

	/* Builtin functions */
	integer s_cmp (char *, char *, ftnlen, ftnlen);

	/* Local variables */
	static integer i__, j, igap;
	static doublereal temp, temp1, temp2;


	/*     %------------------% */
	/*     | Scalar Arguments | */
	/*     %------------------% */


	/*     %-----------------% */
	/*     | Array Arguments | */
	/*     %-----------------% */


	/*     %---------------% */
	/*     | Local Scalars | */
	/*     %---------------% */


	/*     %--------------------% */
	/*     | External Functions | */
	/*     %--------------------% */


	/*     %-----------------------% */
	/*     | Executable Statements | */
	/*     %-----------------------% */

	igap = *n / 2;

	if (s_cmp (which, "LM", (ftnlen) 2, (ftnlen) 2) == 0)
	{

		/*        %------------------------------------------------------% */
		/*        | Sort XREAL,XIMAG into increasing order of magnitude. | */
		/*        %------------------------------------------------------% */

	 L10:
		if (igap == 0)
		{
			goto L9000;
		}

		i__1 = *n - 1;
		for (i__ = igap; i__ <= i__1; ++i__)
		{
			j = i__ - igap;
		 L20:

			if (j < 0)
			{
				goto L30;
			}

			temp1 = dlapy2_ (&xreal[j], &ximag[j]);
			temp2 = dlapy2_ (&xreal[j + igap], &ximag[j + igap]);

			if (temp1 > temp2)
			{
				temp = xreal[j];
				xreal[j] = xreal[j + igap];
				xreal[j + igap] = temp;

				temp = ximag[j];
				ximag[j] = ximag[j + igap];
				ximag[j + igap] = temp;

				if (*apply)
				{
					temp = y[j];
					y[j] = y[j + igap];
					y[j + igap] = temp;
				}
			}
			else
			{
				goto L30;
			}
			j -= igap;
			goto L20;
		 L30:
			;
		}
		igap /= 2;
		goto L10;

	}
	else if (s_cmp (which, "SM", (ftnlen) 2, (ftnlen) 2) == 0)
	{

		/*        %------------------------------------------------------% */
		/*        | Sort XREAL,XIMAG into decreasing order of magnitude. | */
		/*        %------------------------------------------------------% */

	 L40:
		if (igap == 0)
		{
			goto L9000;
		}

		i__1 = *n - 1;
		for (i__ = igap; i__ <= i__1; ++i__)
		{
			j = i__ - igap;
		 L50:

			if (j < 0)
			{
				goto L60;
			}

			temp1 = dlapy2_ (&xreal[j], &ximag[j]);
			temp2 = dlapy2_ (&xreal[j + igap], &ximag[j + igap]);

			if (temp1 < temp2)
			{
				temp = xreal[j];
				xreal[j] = xreal[j + igap];
				xreal[j + igap] = temp;

				temp = ximag[j];
				ximag[j] = ximag[j + igap];
				ximag[j + igap] = temp;

				if (*apply)
				{
					temp = y[j];
					y[j] = y[j + igap];
					y[j + igap] = temp;
				}
			}
			else
			{
				goto L60;
			}
			j -= igap;
			goto L50;
		 L60:
			;
		}
		igap /= 2;
		goto L40;

	}
	else if (s_cmp (which, "LR", (ftnlen) 2, (ftnlen) 2) == 0)
	{

		/*        %------------------------------------------------% */
		/*        | Sort XREAL into increasing order of algebraic. | */
		/*        %------------------------------------------------% */

	 L70:
		if (igap == 0)
		{
			goto L9000;
		}

		i__1 = *n - 1;
		for (i__ = igap; i__ <= i__1; ++i__)
		{
			j = i__ - igap;
		 L80:

			if (j < 0)
			{
				goto L90;
			}

			if (xreal[j] > xreal[j + igap])
			{
				temp = xreal[j];
				xreal[j] = xreal[j + igap];
				xreal[j + igap] = temp;

				temp = ximag[j];
				ximag[j] = ximag[j + igap];
				ximag[j + igap] = temp;

				if (*apply)
				{
					temp = y[j];
					y[j] = y[j + igap];
					y[j + igap] = temp;
				}
			}
			else
			{
				goto L90;
			}
			j -= igap;
			goto L80;
		 L90:
			;
		}
		igap /= 2;
		goto L70;

	}
	else if (s_cmp (which, "SR", (ftnlen) 2, (ftnlen) 2) == 0)
	{

		/*        %------------------------------------------------% */
		/*        | Sort XREAL into decreasing order of algebraic. | */
		/*        %------------------------------------------------% */

	 L100:
		if (igap == 0)
		{
			goto L9000;
		}
		i__1 = *n - 1;
		for (i__ = igap; i__ <= i__1; ++i__)
		{
			j = i__ - igap;
		 L110:

			if (j < 0)
			{
				goto L120;
			}

			if (xreal[j] < xreal[j + igap])
			{
				temp = xreal[j];
				xreal[j] = xreal[j + igap];
				xreal[j + igap] = temp;

				temp = ximag[j];
				ximag[j] = ximag[j + igap];
				ximag[j + igap] = temp;

				if (*apply)
				{
					temp = y[j];
					y[j] = y[j + igap];
					y[j + igap] = temp;
				}
			}
			else
			{
				goto L120;
			}
			j -= igap;
			goto L110;
		 L120:
			;
		}
		igap /= 2;
		goto L100;

	}
	else if (s_cmp (which, "LI", (ftnlen) 2, (ftnlen) 2) == 0)
	{

		/*        %------------------------------------------------% */
		/*        | Sort XIMAG into increasing order of magnitude. | */
		/*        %------------------------------------------------% */

	 L130:
		if (igap == 0)
		{
			goto L9000;
		}
		i__1 = *n - 1;
		for (i__ = igap; i__ <= i__1; ++i__)
		{
			j = i__ - igap;
		 L140:

			if (j < 0)
			{
				goto L150;
			}

			if ((d__1 = ximag[j], abs (d__1)) > (d__2 = ximag[j + igap], abs (d__2)))
			{
				temp = xreal[j];
				xreal[j] = xreal[j + igap];
				xreal[j + igap] = temp;

				temp = ximag[j];
				ximag[j] = ximag[j + igap];
				ximag[j + igap] = temp;

				if (*apply)
				{
					temp = y[j];
					y[j] = y[j + igap];
					y[j + igap] = temp;
				}
			}
			else
			{
				goto L150;
			}
			j -= igap;
			goto L140;
		 L150:
			;
		}
		igap /= 2;
		goto L130;

	}
	else if (s_cmp (which, "SI", (ftnlen) 2, (ftnlen) 2) == 0)
	{

		/*        %------------------------------------------------% */
		/*        | Sort XIMAG into decreasing order of magnitude. | */
		/*        %------------------------------------------------% */

	 L160:
		if (igap == 0)
		{
			goto L9000;
		}
		i__1 = *n - 1;
		for (i__ = igap; i__ <= i__1; ++i__)
		{
			j = i__ - igap;
		 L170:

			if (j < 0)
			{
				goto L180;
			}

			if ((d__1 = ximag[j], abs (d__1)) < (d__2 = ximag[j + igap], abs (d__2)))
			{
				temp = xreal[j];
				xreal[j] = xreal[j + igap];
				xreal[j + igap] = temp;

				temp = ximag[j];
				ximag[j] = ximag[j + igap];
				ximag[j + igap] = temp;

				if (*apply)
				{
					temp = y[j];
					y[j] = y[j + igap];
					y[j + igap] = temp;
				}
			}
			else
			{
				goto L180;
			}
			j -= igap;
			goto L170;
		 L180:
			;
		}
		igap /= 2;
		goto L160;
	}

 L9000:
	return 0;

	/*     %---------------% */
	/*     | End of dsortc | */
	/*     %---------------% */

}	/* dsortc_ */


/*     %---------------------------------------------% */
/*     | Initialize statistic and timing information | */
/*     | for nonsymmetric Arnoldi code.              | */
/*     %---------------------------------------------% */

/* \Author */
/*     Danny Sorensen               Phuong Vu */
/*     Richard Lehoucq              CRPC / Rice University */
/*     Dept. of Computational &     Houston, Texas */
/*     Applied Mathematics */
/*     Rice University */
/*     Houston, Texas */

/* \SCCS Information: @(#) */
/* FILE: statn.F   SID: 2.4   DATE OF SID: 4/20/96   RELEASE: 2 */

static int
dstatn_ (void)
{

	/*     %--------------------------------% */
	/*     | See stat.doc for documentation | */
	/*     %--------------------------------% */


	/*     %-----------------------% */
	/*     | Executable Statements | */
	/*     %-----------------------% */

	/*     %--------------------------------% */
	/*     | See stat.doc for documentation | */
	/*     %--------------------------------% */

	/* \SCCS Information: @(#) */
	/* FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 */


	timing_1.nopx = 0;
	timing_1.nbx = 0;
	timing_1.nrorth = 0;
	timing_1.nitref = 0;
	timing_1.nrstrt = 0;

	timing_1.tnaupd = 0.f;
	timing_1.tnaup2 = 0.f;
	timing_1.tnaitr = 0.f;
	timing_1.tneigh = 0.f;
	timing_1.tngets = 0.f;
	timing_1.tnapps = 0.f;
	timing_1.tnconv = 0.f;
	timing_1.titref = 0.f;
	timing_1.tgetv0 = 0.f;
	timing_1.trvec = 0.f;

	/*     %----------------------------------------------------% */
	/*     | User time including reverse communication overhead | */
	/*     %----------------------------------------------------% */

	timing_1.tmvopx = 0.f;
	timing_1.tmvbx = 0.f;

	return 0;


	/*     %---------------% */
	/*     | End of dstatn | */
	/*     %---------------% */

}	/* dstatn_ */

/* ----------------------------------------------------------------------- */
/*  Routine:    DVOUT */

/*  Purpose:    Real vector output routine. */

/*  Usage:      CALL DVOUT (LOUT, N, SX, IDIGIT, IFMT) */

/*  Arguments */
/*     N      - Length of array SX.  (Input) */
/*     SX     - Real array to be printed.  (Input) */
/*     IFMT   - Format to be used in printing array SX.  (Input) */
/*     IDIGIT - Print up to IABS(IDIGIT) decimal digits per number.  (In) */
/*              If IDIGIT .LT. 0, printing is done with 72 columns. */
/*              If IDIGIT .GT. 0, printing is done with 132 columns. */

/* ----------------------------------------------------------------------- */

static int
dvout_ (integer * lout, integer * n, doublereal * sx, integer * idigit, char *ifmt, ftnlen ifmt_len)
{
	/* Format strings */
	static char fmt_9999[] = "(/1x,a,/1x,a)";
	static char fmt_9998[] = "(1x,i4,\002 - \002,i4,\002:\002,1p,10d12.3)";
	static char fmt_9997[] = "(1x,i4,\002 - \002,i4,\002:\002,1x,1p,8d14.5)";
	static char fmt_9996[] = "(1x,i4,\002 - \002,i4,\002:\002,1x,1p,6d18.9)";
	static char fmt_9995[] = "(1x,i4,\002 - \002,i4,\002:\002,1x,1p,5d24.13)";
	static char fmt_9994[] = "(1x,\002 \002)";

	/* System generated locals */
	integer i__1, i__2, i__3;

	/* Builtin functions */
	integer i_len (char *, ftnlen), s_wsfe (cilist *), do_fio (integer *, char *, ftnlen), e_wsfe (void);

	/* Local variables */
	static integer i__, k1, k2, lll;
	static char line[80];
	static integer ndigit;

	/* Fortran I/O blocks */
	static cilist io___1067 = { 0, 0, 0, fmt_9999, 0 };
	static cilist io___1071 = { 0, 0, 0, fmt_9998, 0 };
	static cilist io___1072 = { 0, 0, 0, fmt_9997, 0 };
	static cilist io___1073 = { 0, 0, 0, fmt_9996, 0 };
	static cilist io___1074 = { 0, 0, 0, fmt_9995, 0 };
	static cilist io___1075 = { 0, 0, 0, fmt_9998, 0 };
	static cilist io___1076 = { 0, 0, 0, fmt_9997, 0 };
	static cilist io___1077 = { 0, 0, 0, fmt_9996, 0 };
	static cilist io___1078 = { 0, 0, 0, fmt_9995, 0 };
	static cilist io___1079 = { 0, 0, 0, fmt_9994, 0 };


	/*     ... */
	/*     ... SPECIFICATIONS FOR ARGUMENTS */
	/*     ... */
	/*     ... SPECIFICATIONS FOR LOCAL VARIABLES */
	/*     .. Scalar Arguments .. */
	/*     .. */
	/*     .. Array Arguments .. */
	/*     .. */
	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Intrinsic Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */
	/*     ... */
	/*     ... FIRST EXECUTABLE STATEMENT */


	/* Parameter adjustments */
	--sx;

	/* Function Body */
	/* Computing MIN */
	i__1 = i_len (ifmt, ifmt_len);
	lll = min (i__1, 80);
	i__1 = lll;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		*(unsigned char *) &line[i__ - 1] = '-';
		/* L10: */
	}

	for (i__ = lll + 1; i__ <= 80; ++i__)
	{
		*(unsigned char *) &line[i__ - 1] = ' ';
		/* L20: */
	}

	io___1067.ciunit = *lout;
	s_wsfe (&io___1067);
	do_fio (&c__1, ifmt, ifmt_len);
	do_fio (&c__1, line, lll);
	e_wsfe ();

	if (*n <= 0)
	{
		return 0;
	}
	ndigit = *idigit;
	if (*idigit == 0)
	{
		ndigit = 4;
	}

	/* ======================================================================= */
	/*             CODE FOR OUTPUT USING 72 COLUMNS FORMAT */
	/* ======================================================================= */

	if (*idigit < 0)
	{
		ndigit = -(*idigit);
		if (ndigit <= 4)
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 5)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 4;
				k2 = min (i__2, i__3);
				io___1071.ciunit = *lout;
				s_wsfe (&io___1071);
				do_fio (&c__1, (char *) &k1, (ftnlen) sizeof (integer));
				do_fio (&c__1, (char *) &k2, (ftnlen) sizeof (integer));
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__1, (char *) &sx[i__], (ftnlen) sizeof (doublereal));
				}
				e_wsfe ();
				/* L30: */
			}
		}
		else if (ndigit <= 6)
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 4)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 3;
				k2 = min (i__2, i__3);
				io___1072.ciunit = *lout;
				s_wsfe (&io___1072);
				do_fio (&c__1, (char *) &k1, (ftnlen) sizeof (integer));
				do_fio (&c__1, (char *) &k2, (ftnlen) sizeof (integer));
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__1, (char *) &sx[i__], (ftnlen) sizeof (doublereal));
				}
				e_wsfe ();
				/* L40: */
			}
		}
		else if (ndigit <= 10)
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 3)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 2;
				k2 = min (i__2, i__3);
				io___1073.ciunit = *lout;
				s_wsfe (&io___1073);
				do_fio (&c__1, (char *) &k1, (ftnlen) sizeof (integer));
				do_fio (&c__1, (char *) &k2, (ftnlen) sizeof (integer));
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__1, (char *) &sx[i__], (ftnlen) sizeof (doublereal));
				}
				e_wsfe ();
				/* L50: */
			}
		}
		else
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 2)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 1;
				k2 = min (i__2, i__3);
				io___1074.ciunit = *lout;
				s_wsfe (&io___1074);
				do_fio (&c__1, (char *) &k1, (ftnlen) sizeof (integer));
				do_fio (&c__1, (char *) &k2, (ftnlen) sizeof (integer));
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__1, (char *) &sx[i__], (ftnlen) sizeof (doublereal));
				}
				e_wsfe ();
				/* L60: */
			}
		}

		/* ======================================================================= */
		/*             CODE FOR OUTPUT USING 132 COLUMNS FORMAT */
		/* ======================================================================= */

	}
	else
	{
		if (ndigit <= 4)
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 10)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 9;
				k2 = min (i__2, i__3);
				io___1075.ciunit = *lout;
				s_wsfe (&io___1075);
				do_fio (&c__1, (char *) &k1, (ftnlen) sizeof (integer));
				do_fio (&c__1, (char *) &k2, (ftnlen) sizeof (integer));
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__1, (char *) &sx[i__], (ftnlen) sizeof (doublereal));
				}
				e_wsfe ();
				/* L70: */
			}
		}
		else if (ndigit <= 6)
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 8)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 7;
				k2 = min (i__2, i__3);
				io___1076.ciunit = *lout;
				s_wsfe (&io___1076);
				do_fio (&c__1, (char *) &k1, (ftnlen) sizeof (integer));
				do_fio (&c__1, (char *) &k2, (ftnlen) sizeof (integer));
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__1, (char *) &sx[i__], (ftnlen) sizeof (doublereal));
				}
				e_wsfe ();
				/* L80: */
			}
		}
		else if (ndigit <= 10)
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 6)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 5;
				k2 = min (i__2, i__3);
				io___1077.ciunit = *lout;
				s_wsfe (&io___1077);
				do_fio (&c__1, (char *) &k1, (ftnlen) sizeof (integer));
				do_fio (&c__1, (char *) &k2, (ftnlen) sizeof (integer));
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__1, (char *) &sx[i__], (ftnlen) sizeof (doublereal));
				}
				e_wsfe ();
				/* L90: */
			}
		}
		else
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 5)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 4;
				k2 = min (i__2, i__3);
				io___1078.ciunit = *lout;
				s_wsfe (&io___1078);
				do_fio (&c__1, (char *) &k1, (ftnlen) sizeof (integer));
				do_fio (&c__1, (char *) &k2, (ftnlen) sizeof (integer));
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__1, (char *) &sx[i__], (ftnlen) sizeof (doublereal));
				}
				e_wsfe ();
				/* L100: */
			}
		}
	}
	io___1079.ciunit = *lout;
	s_wsfe (&io___1079);
	e_wsfe ();
	return 0;
}	/* dvout_ */

/* ----------------------------------------------------------------------- */
/*  Routine:    IVOUT */

/*  Purpose:    Integer vector output routine. */

/*  Usage:      CALL IVOUT (LOUT, N, IX, IDIGIT, IFMT) */

/*  Arguments */
/*     N      - Length of array IX. (Input) */
/*     IX     - Integer array to be printed. (Input) */
/*     IFMT   - Format to be used in printing array IX. (Input) */
/*     IDIGIT - Print up to ABS(IDIGIT) decimal digits / number. (Input) */
/*              If IDIGIT .LT. 0, printing is done with 72 columns. */
/*              If IDIGIT .GT. 0, printing is done with 132 columns. */

/* ----------------------------------------------------------------------- */

static int
ivout_ (integer * lout, integer * n, integer * ix, integer * idigit, char *ifmt, ftnlen ifmt_len)
{
	/* Format strings */
	static char fmt_2000[] = "(/1x,a/1x,a)";
	static char fmt_1000[] = "(1x,i4,\002 - \002,i4,\002:\002,20(1x,i5))";
	static char fmt_1001[] = "(1x,i4,\002 - \002,i4,\002:\002,15(1x,i7))";
	static char fmt_1002[] = "(1x,i4,\002 - \002,i4,\002:\002,10(1x,i11))";
	static char fmt_1003[] = "(1x,i4,\002 - \002,i4,\002:\002,7(1x,i15))";
	static char fmt_1004[] = "(1x,\002 \002)";

	/* System generated locals */
	integer i__1, i__2, i__3;

	/* Builtin functions */
	integer i_len (char *, ftnlen), s_wsfe (cilist *), do_fio (integer *, char *, ftnlen), e_wsfe (void);

	/* Local variables */
	static integer i__, k1, k2, lll;
	static char line[80];
	static integer ndigit;

	/* Fortran I/O blocks */
	static cilist io___1083 = { 0, 0, 0, fmt_2000, 0 };
	static cilist io___1087 = { 0, 0, 0, fmt_1000, 0 };
	static cilist io___1088 = { 0, 0, 0, fmt_1001, 0 };
	static cilist io___1089 = { 0, 0, 0, fmt_1002, 0 };
	static cilist io___1090 = { 0, 0, 0, fmt_1003, 0 };
	static cilist io___1091 = { 0, 0, 0, fmt_1000, 0 };
	static cilist io___1092 = { 0, 0, 0, fmt_1001, 0 };
	static cilist io___1093 = { 0, 0, 0, fmt_1002, 0 };
	static cilist io___1094 = { 0, 0, 0, fmt_1003, 0 };
	static cilist io___1095 = { 0, 0, 0, fmt_1004, 0 };


	/*     ... */
	/*     ... SPECIFICATIONS FOR ARGUMENTS */
	/*     ... */
	/*     ... SPECIFICATIONS FOR LOCAL VARIABLES */
	/*     ... */
	/*     ... SPECIFICATIONS INTRINSICS */


	/* Parameter adjustments */
	--ix;

	/* Function Body */
	/* Computing MIN */
	i__1 = i_len (ifmt, ifmt_len);
	lll = min (i__1, 80);
	i__1 = lll;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		*(unsigned char *) &line[i__ - 1] = '-';
		/* L1: */
	}

	for (i__ = lll + 1; i__ <= 80; ++i__)
	{
		*(unsigned char *) &line[i__ - 1] = ' ';
		/* L2: */
	}

	io___1083.ciunit = *lout;
	s_wsfe (&io___1083);
	do_fio (&c__1, ifmt, ifmt_len);
	do_fio (&c__1, line, lll);
	e_wsfe ();

	if (*n <= 0)
	{
		return 0;
	}
	ndigit = *idigit;
	if (*idigit == 0)
	{
		ndigit = 4;
	}

	/* ======================================================================= */
	/*             CODE FOR OUTPUT USING 72 COLUMNS FORMAT */
	/* ======================================================================= */

	if (*idigit < 0)
	{

		ndigit = -(*idigit);
		if (ndigit <= 4)
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 10)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 9;
				k2 = min (i__2, i__3);
				io___1087.ciunit = *lout;
				s_wsfe (&io___1087);
				do_fio (&c__1, (char *) &k1, (ftnlen) sizeof (integer));
				do_fio (&c__1, (char *) &k2, (ftnlen) sizeof (integer));
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__1, (char *) &ix[i__], (ftnlen) sizeof (integer));
				}
				e_wsfe ();
				/* L10: */
			}

		}
		else if (ndigit <= 6)
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 7)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 6;
				k2 = min (i__2, i__3);
				io___1088.ciunit = *lout;
				s_wsfe (&io___1088);
				do_fio (&c__1, (char *) &k1, (ftnlen) sizeof (integer));
				do_fio (&c__1, (char *) &k2, (ftnlen) sizeof (integer));
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__1, (char *) &ix[i__], (ftnlen) sizeof (integer));
				}
				e_wsfe ();
				/* L30: */
			}

		}
		else if (ndigit <= 10)
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 5)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 4;
				k2 = min (i__2, i__3);
				io___1089.ciunit = *lout;
				s_wsfe (&io___1089);
				do_fio (&c__1, (char *) &k1, (ftnlen) sizeof (integer));
				do_fio (&c__1, (char *) &k2, (ftnlen) sizeof (integer));
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__1, (char *) &ix[i__], (ftnlen) sizeof (integer));
				}
				e_wsfe ();
				/* L50: */
			}

		}
		else
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 3)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 2;
				k2 = min (i__2, i__3);
				io___1090.ciunit = *lout;
				s_wsfe (&io___1090);
				do_fio (&c__1, (char *) &k1, (ftnlen) sizeof (integer));
				do_fio (&c__1, (char *) &k2, (ftnlen) sizeof (integer));
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__1, (char *) &ix[i__], (ftnlen) sizeof (integer));
				}
				e_wsfe ();
				/* L70: */
			}
		}

		/* ======================================================================= */
		/*             CODE FOR OUTPUT USING 132 COLUMNS FORMAT */
		/* ======================================================================= */

	}
	else
	{

		if (ndigit <= 4)
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 20)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 19;
				k2 = min (i__2, i__3);
				io___1091.ciunit = *lout;
				s_wsfe (&io___1091);
				do_fio (&c__1, (char *) &k1, (ftnlen) sizeof (integer));
				do_fio (&c__1, (char *) &k2, (ftnlen) sizeof (integer));
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__1, (char *) &ix[i__], (ftnlen) sizeof (integer));
				}
				e_wsfe ();
				/* L90: */
			}

		}
		else if (ndigit <= 6)
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 15)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 14;
				k2 = min (i__2, i__3);
				io___1092.ciunit = *lout;
				s_wsfe (&io___1092);
				do_fio (&c__1, (char *) &k1, (ftnlen) sizeof (integer));
				do_fio (&c__1, (char *) &k2, (ftnlen) sizeof (integer));
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__1, (char *) &ix[i__], (ftnlen) sizeof (integer));
				}
				e_wsfe ();
				/* L110: */
			}

		}
		else if (ndigit <= 10)
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 10)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 9;
				k2 = min (i__2, i__3);
				io___1093.ciunit = *lout;
				s_wsfe (&io___1093);
				do_fio (&c__1, (char *) &k1, (ftnlen) sizeof (integer));
				do_fio (&c__1, (char *) &k2, (ftnlen) sizeof (integer));
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__1, (char *) &ix[i__], (ftnlen) sizeof (integer));
				}
				e_wsfe ();
				/* L130: */
			}

		}
		else
		{
			i__1 = *n;
			for (k1 = 1; k1 <= i__1; k1 += 7)
			{
				/* Computing MIN */
				i__2 = *n, i__3 = k1 + 6;
				k2 = min (i__2, i__3);
				io___1094.ciunit = *lout;
				s_wsfe (&io___1094);
				do_fio (&c__1, (char *) &k1, (ftnlen) sizeof (integer));
				do_fio (&c__1, (char *) &k2, (ftnlen) sizeof (integer));
				i__2 = k2;
				for (i__ = k1; i__ <= i__2; ++i__)
				{
					do_fio (&c__1, (char *) &ix[i__], (ftnlen) sizeof (integer));
				}
				e_wsfe ();
				/* L150: */
			}
		}
	}
	io___1095.ciunit = *lout;
	s_wsfe (&io___1095);
	e_wsfe ();


	return 0;
}	/* ivout_ */

static int
second_ (real * t)
{
	static real t1;
	static real tarray[2];



	/*  -- LAPACK auxiliary routine (preliminary version) -- */
	/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
	/*     Courant Institute, Argonne National Lab, and Rice University */
	/*     July 26, 1991 */

	/*  Purpose */
	/*  ======= */

	/*  SECOND returns the user time for a process in seconds. */
	/*  This version gets the time from the system function ETIME. */

	/*     .. Local Scalars .. */
	/*     .. */
	/*     .. Local Arrays .. */
	/*     .. */
	/*     .. External Functions .. */
	/*     .. */
	/*     .. Executable Statements .. */

	t1 = etime_ (tarray);
	*t = tarray[0];
	return 0;

	/*     End of SECOND */

}	/* second_ */
