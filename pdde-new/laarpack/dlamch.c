#include "cblaswrap.h"

double dlamch_(char *cmach, int cmach_len);

static int dlamc5_(int * beta, int * p, int * emin, bool * ieee, int * emax, double * rmax);
static int dlamc4_(int * emin, double * start, int * base);
static double dlamc3_(double * a, double * b);
static int dlamc1_(int * beta, int * t, bool * rnd, bool * ieee1);
static int dlamc2_(int * beta, int * t, bool * rnd, double * eps, int * emin, double * rmin, int * emax, double * rmax);

/* Table of constant values */

static int c__1 = 1;
static double c_b507 = 0.;

double
dlamch_(char *cmach, int cmach_len)
{
  /* Initialized data */

  static bool first = true;

  /* System generated locals */
  int i__1;
  double ret_val;

  /* Local variables */
  static double t;
  static int it;
  static double rnd, eps, base;
  static int beta;
  static double emin, prec, emax;
  static int imin, imax;
  static bool lrnd;
  static double rmin, rmax, rmach;
  static double small, sfmin;


  /*  -- LAPACK auxiliary routine (version 3.0) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLAMCH determines double precision machine parameters. */

  /*  Arguments */
  /*  ========= */

  /*  CMACH   (input) CHARACTER*1 */
  /*          Specifies the value to be returned by DLAMCH: */
  /*          = 'E' or 'e',   DLAMCH := eps */
  /*          = 'S' or 's ,   DLAMCH := sfmin */
  /*          = 'B' or 'b',   DLAMCH := base */
  /*          = 'P' or 'p',   DLAMCH := eps*base */
  /*          = 'N' or 'n',   DLAMCH := t */
  /*          = 'R' or 'r',   DLAMCH := rnd */
  /*          = 'M' or 'm',   DLAMCH := emin */
  /*          = 'U' or 'u',   DLAMCH := rmin */
  /*          = 'L' or 'l',   DLAMCH := emax */
  /*          = 'O' or 'o',   DLAMCH := rmax */

  /*          where */

  /*          eps   = relative machine precision */
  /*          sfmin = safe minimum, such that 1/sfmin does not overflow */
  /*          base  = base of the machine */
  /*          prec  = eps*base */
  /*          t     = number of (base) digits in the mantissa */
  /*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise */
  /*          emin  = minimum exponent before (gradual) underflow */
  /*          rmin  = underflow threshold - base**(emin-1) */
  /*          emax  = largest exponent before overflow */
  /*          rmax  = overflow threshold  - (base**emax)*(1-eps) */

  /* ===================================================================== */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. External Subroutines .. */
  /*     .. */
  /*     .. Save statement .. */
  /*     .. */
  /*     .. Data statements .. */
  /*     .. */
  /*     .. Executable Statements .. */

  if (first)
  {
    first = false;
    dlamc2_(&beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax);
    base = (double) beta;
    t = (double) it;
    if (lrnd)
    {
      rnd = 1.;
      i__1 = 1 - it;
      eps = pow_di(&base, &i__1) / 2;
    }
    else
    {
      rnd = 0.;
      i__1 = 1 - it;
      eps = pow_di(&base, &i__1);
    }
    prec = eps * base;
    emin = (double) imin;
    emax = (double) imax;
    sfmin = rmin;
    small = 1. / rmax;
    if (small >= sfmin)
    {

      /*           Use SMALL plus a bit, to avoid the possibility of rounding */
      /*           causing overflow when computing  1/sfmin. */

      sfmin = small * (eps + 1.);
    }
  }

  if (lsame_(cmach, "E", (int) 1, (int) 1))
  {
    rmach = eps;
  }
  else if (lsame_(cmach, "S", (int) 1, (int) 1))
  {
    rmach = sfmin;
  }
  else if (lsame_(cmach, "B", (int) 1, (int) 1))
  {
    rmach = base;
  }
  else if (lsame_(cmach, "P", (int) 1, (int) 1))
  {
    rmach = prec;
  }
  else if (lsame_(cmach, "N", (int) 1, (int) 1))
  {
    rmach = t;
  }
  else if (lsame_(cmach, "R", (int) 1, (int) 1))
  {
    rmach = rnd;
  }
  else if (lsame_(cmach, "M", (int) 1, (int) 1))
  {
    rmach = emin;
  }
  else if (lsame_(cmach, "U", (int) 1, (int) 1))
  {
    rmach = rmin;
  }
  else if (lsame_(cmach, "L", (int) 1, (int) 1))
  {
    rmach = emax;
  }
  else if (lsame_(cmach, "O", (int) 1, (int) 1))
  {
    rmach = rmax;
  }

  ret_val = rmach;
  return ret_val;

  /*     End of DLAMCH */

} /* dlamch_ */


/* *********************************************************************** */

static int
dlamc1_(int * beta, int * t, bool * rnd, bool * ieee1)
{
  /* Initialized data */

  static bool first = true;

  /* System generated locals */
  double d__1, d__2;

  /* Local variables */
  static double a, b, c__, f, t1, t2;
  static int lt;
  static double one, qtr;
  static bool lrnd;
  static int lbeta;
  static double savec;
  static bool lieee1;


  /*  -- LAPACK auxiliary routine (version 3.0) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLAMC1 determines the machine parameters given by BETA, T, RND, and */
  /*  IEEE1. */

  /*  Arguments */
  /*  ========= */

  /*  BETA    (output) INTEGER */
  /*          The base of the machine. */

  /*  T       (output) INTEGER */
  /*          The number of ( BETA ) digits in the mantissa. */

  /*  RND     (output) LOGICAL */
  /*          Specifies whether proper rounding  ( RND = .TRUE. )  or */
  /*          chopping  ( RND = .FALSE. )  occurs in addition. This may not */
  /*          be a reliable guide to the way in which the machine performs */
  /*          its arithmetic. */

  /*  IEEE1   (output) LOGICAL */
  /*          Specifies whether rounding appears to be done in the IEEE */
  /*          'round to nearest' style. */

  /*  Further Details */
  /*  =============== */

  /*  The routine is based on the routine  ENVRON  by Malcolm and */
  /*  incorporates suggestions by Gentleman and Marovich. See */

  /*     Malcolm M. A. (1972) Algorithms to reveal properties of */
  /*        floating-point arithmetic. Comms. of the ACM, 15, 949-951. */

  /*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms */
  /*        that reveal properties of floating point arithmetic units. */
  /*        Comms. of the ACM, 17, 276-277. */

  /* ===================================================================== */

  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Save statement .. */
  /*     .. */
  /*     .. Data statements .. */
  /*     .. */
  /*     .. Executable Statements .. */

  if (first)
  {
    first = false;
    one = 1.;

    /*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA, */
    /*        IEEE1, T and RND. */

    /*        Throughout this routine  we use the function  DLAMC3  to ensure */
    /*        that relevant values are  stored and not held in registers,  or */
    /*        are not affected by optimizers. */

    /*        Compute  a = 2.0**m  with the  smallest positive int m such */
    /*        that */

    /*           fl( a + 1.0 ) = a. */

    a = 1.;
    c__ = 1.;

    /* +       WHILE( C.EQ.ONE )LOOP */
L10:
    if (c__ == one)
    {
      a *= 2;
      c__ = dlamc3_(&a, &one);
      d__1 = -a;
      c__ = dlamc3_(&c__, &d__1);
      goto L10;
    }
    /* +       END WHILE */

    /*        Now compute  b = 2.0**m  with the smallest positive int m */
    /*        such that */

    /*           fl( a + b ) .gt. a. */

    b = 1.;
    c__ = dlamc3_(&a, &b);

    /* +       WHILE( C.EQ.A )LOOP */
L20:
    if (c__ == a)
    {
      b *= 2;
      c__ = dlamc3_(&a, &b);
      goto L20;
    }
    /* +       END WHILE */

    /*        Now compute the base.  a and c  are neighbouring floating point */
    /*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so */
    /*        their difference is beta. Adding 0.25 to c is to ensure that it */
    /*        is truncated to beta and not ( beta - 1 ). */

    qtr = one / 4;
    savec = c__;
    d__1 = -a;
    c__ = dlamc3_(&c__, &d__1);
    lbeta = (int)(c__ + qtr);

    /*        Now determine whether rounding or chopping occurs,  by adding a */
    /*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a. */

    b = (double) lbeta;
    d__1 = b / 2;
    d__2 = -b / 100;
    f = dlamc3_(&d__1, &d__2);
    c__ = dlamc3_(&f, &a);
    if (c__ == a)
    {
      lrnd = true;
    }
    else
    {
      lrnd = false;
    }
    d__1 = b / 2;
    d__2 = b / 100;
    f = dlamc3_(&d__1, &d__2);
    c__ = dlamc3_(&f, &a);
    if (lrnd && c__ == a)
    {
      lrnd = false;
    }

    /*        Try and decide whether rounding is done in the  IEEE  'round to */
    /*        nearest' style. B/2 is half a unit in the last place of the two */
    /*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit */
    /*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change */
    /*        A, but adding B/2 to SAVEC should change SAVEC. */

    d__1 = b / 2;
    t1 = dlamc3_(&d__1, &a);
    d__1 = b / 2;
    t2 = dlamc3_(&d__1, &savec);
    lieee1 = t1 == a && t2 > savec && lrnd;

    /*        Now find  the  mantissa, t.  It should  be the  int part of */
    /*        log to the base beta of a,  however it is safer to determine  t */
    /*        by powering.  So we find t as the smallest positive int for */
    /*        which */

    /*           fl( beta**t + 1.0 ) = 1.0. */

    lt = 0;
    a = 1.;
    c__ = 1.;

    /* +       WHILE( C.EQ.ONE )LOOP */
L30:
    if (c__ == one)
    {
      ++lt;
      a *= lbeta;
      c__ = dlamc3_(&a, &one);
      d__1 = -a;
      c__ = dlamc3_(&c__, &d__1);
      goto L30;
    }
    /* +       END WHILE */

  }

  *beta = lbeta;
  *t = lt;
  *rnd = lrnd;
  *ieee1 = lieee1;
  return 0;

  /*     End of DLAMC1 */

} /* dlamc1_ */


/* *********************************************************************** */

static int
dlamc2_(int * beta, int * t, bool * rnd, double * eps, int * emin, double * rmin, int * emax, double * rmax)
{
  /* Initialized data */

  static bool first = true;
  static bool iwarn = false;

  /* System generated locals */
  int i__1;
  double d__1, d__2, d__3, d__4, d__5;

  /* Local variables */
  static double a, b, c__;
  static int i__, lt;
  static double one, two;
  static bool ieee;
  static double half;
  static bool lrnd;
  static double leps, zero;
  static int lbeta;
  static double rbase;
  static int lemin, lemax, gnmin;
  static double small;
  static int gpmin;
  static double third, lrmin, lrmax, sixth;
  static bool lieee1;
  static int ngnmin, ngpmin;

  /*  -- LAPACK auxiliary routine (version 3.0) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLAMC2 determines the machine parameters specified in its argument */
  /*  list. */

  /*  Arguments */
  /*  ========= */

  /*  BETA    (output) INTEGER */
  /*          The base of the machine. */

  /*  T       (output) INTEGER */
  /*          The number of ( BETA ) digits in the mantissa. */

  /*  RND     (output) LOGICAL */
  /*          Specifies whether proper rounding  ( RND = .TRUE. )  or */
  /*          chopping  ( RND = .FALSE. )  occurs in addition. This may not */
  /*          be a reliable guide to the way in which the machine performs */
  /*          its arithmetic. */

  /*  EPS     (output) DOUBLE PRECISION */
  /*          The smallest positive number such that */

  /*             fl( 1.0 - EPS ) .LT. 1.0, */

  /*          where fl denotes the computed value. */

  /*  EMIN    (output) INTEGER */
  /*          The minimum exponent before (gradual) underflow occurs. */

  /*  RMIN    (output) DOUBLE PRECISION */
  /*          The smallest normalized number for the machine, given by */
  /*          BASE**( EMIN - 1 ), where  BASE  is the floating point value */
  /*          of BETA. */

  /*  EMAX    (output) INTEGER */
  /*          The maximum exponent before overflow occurs. */

  /*  RMAX    (output) DOUBLE PRECISION */
  /*          The largest positive number for the machine, given by */
  /*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point */
  /*          value of BETA. */

  /*  Further Details */
  /*  =============== */

  /*  The computation of  EPS  is based on a routine PARANOIA by */
  /*  W. Kahan of the University of California at Berkeley. */

  /* ===================================================================== */

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
  /*     .. Data statements .. */
  /*     .. */
  /*     .. Executable Statements .. */

  if (first)
  {
    first = false;
    zero = 0.;
    one = 1.;
    two = 2.;

    /*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of */
    /*        BETA, T, RND, EPS, EMIN and RMIN. */

    /*        Throughout this routine  we use the function  DLAMC3  to ensure */
    /*        that relevant values are stored  and not held in registers,  or */
    /*        are not affected by optimizers. */

    /*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1. */

    dlamc1_(&lbeta, &lt, &lrnd, &lieee1);

    /*        Start to find EPS. */

    b = (double) lbeta;
    i__1 = -lt;
    a = pow_di(&b, &i__1);
    leps = a;

    /*        Try some tricks to see whether or not this is the correct  EPS. */

    b = two / 3;
    half = one / 2;
    d__1 = -half;
    sixth = dlamc3_(&b, &d__1);
    third = dlamc3_(&sixth, &sixth);
    d__1 = -half;
    b = dlamc3_(&third, &d__1);
    b = dlamc3_(&b, &sixth);
    b = abs(b);
    if (b < leps)
    {
      b = leps;
    }

    leps = 1.;

    /* +       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP */
L10:
    if (leps > b && b > zero)
    {
      leps = b;
      d__1 = half * leps;
      /* Computing 5th power */
      d__3 = two, d__4 = d__3, d__3 *= d__3;
      /* Computing 2nd power */
      d__5 = leps;
      d__2 = d__4 * (d__3 * d__3) * (d__5 * d__5);
      c__ = dlamc3_(&d__1, &d__2);
      d__1 = -c__;
      c__ = dlamc3_(&half, &d__1);
      b = dlamc3_(&half, &c__);
      d__1 = -b;
      c__ = dlamc3_(&half, &d__1);
      b = dlamc3_(&half, &c__);
      goto L10;
    }
    /* +       END WHILE */

    if (a < leps)
    {
      leps = a;
    }

    /*        Computation of EPS complete. */

    /*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)). */
    /*        Keep dividing  A by BETA until (gradual) underflow occurs. This */
    /*        is detected when we cannot recover the previous A. */

    rbase = one / lbeta;
    small = one;
    for (i__ = 1; i__ <= 3; ++i__)
    {
      d__1 = small * rbase;
      small = dlamc3_(&d__1, &zero);
      /* L20: */
    }
    a = dlamc3_(&one, &small);
    dlamc4_(&ngpmin, &one, &lbeta);
    d__1 = -one;
    dlamc4_(&ngnmin, &d__1, &lbeta);
    dlamc4_(&gpmin, &a, &lbeta);
    d__1 = -a;
    dlamc4_(&gnmin, &d__1, &lbeta);
    ieee = false;

    if (ngpmin == ngnmin && gpmin == gnmin)
    {
      if (ngpmin == gpmin)
      {
        lemin = ngpmin;
        /*            ( Non twos-complement machines, no gradual underflow; */
        /*              e.g.,  VAX ) */
      }
      else if (gpmin - ngpmin == 3)
      {
        lemin = ngpmin - 1 + lt;
        ieee = true;
        /*            ( Non twos-complement machines, with gradual underflow; */
        /*              e.g., IEEE standard followers ) */
      }
      else
      {
        lemin = min(ngpmin, gpmin);
        /*            ( A guess; no known machine ) */
        iwarn = true;
      }

    }
    else if (ngpmin == gpmin && ngnmin == gnmin)
    {
      if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1)
      {
        lemin = max(ngpmin, ngnmin);
        /*            ( Twos-complement machines, no gradual underflow; */
        /*              e.g., CYBER 205 ) */
      }
      else
      {
        lemin = min(ngpmin, ngnmin);
        /*            ( A guess; no known machine ) */
        iwarn = true;
      }

    }
    else if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1 && gpmin == gnmin)
    {
      if (gpmin - min(ngpmin, ngnmin) == 3)
      {
        lemin = max(ngpmin, ngnmin) - 1 + lt;
        /*            ( Twos-complement machines with gradual underflow; */
        /*              no known machine ) */
      }
      else
      {
        lemin = min(ngpmin, ngnmin);
        /*            ( A guess; no known machine ) */
        iwarn = true;
      }

    }
    else
    {
      /* Computing MIN */
      i__1 = min(ngpmin, ngnmin), i__1 = min(i__1, gpmin);
      lemin = min(i__1, gnmin);
      /*         ( A guess; no known machine ) */
      iwarn = true;
    }
    /* ** */
    /* Comment out this if block if EMIN is ok */
    if (iwarn)
    {
//    first = true;
//    s_wsfe (&io___381);
//    do_fio (&c__1, (char *) &lemin, (int) sizeof (int));
//    e_wsfe ();
    }
    /* ** */

    /*        Assume IEEE arithmetic if we found denormalised  numbers above, */
    /*        or if arithmetic seems to round in the  IEEE style,  determined */
    /*        in routine DLAMC1. A true IEEE machine should have both  things */
    /*        true; however, faulty machines may have one or the other. */

    ieee = ieee || lieee1;

    /*        Compute  RMIN by successive division by  BETA. We could compute */
    /*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during */
    /*        this computation. */

    lrmin = 1.;
    i__1 = 1 - lemin;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      d__1 = lrmin * rbase;
      lrmin = dlamc3_(&d__1, &zero);
      /* L30: */
    }

    /*        Finally, call DLAMC5 to compute EMAX and RMAX. */

    dlamc5_(&lbeta, &lt, &lemin, &ieee, &lemax, &lrmax);
  }

  *beta = lbeta;
  *t = lt;
  *rnd = lrnd;
  *eps = leps;
  *emin = lemin;
  *rmin = lrmin;
  *emax = lemax;
  *rmax = lrmax;

  return 0;


  /*     End of DLAMC2 */

} /* dlamc2_ */


/* *********************************************************************** */

static double
dlamc3_(double * a, double * b)
{
  /* System generated locals */
  double ret_val;


  /*  -- LAPACK auxiliary routine (version 3.0) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing */
  /*  the addition of  A  and  B ,  for use in situations where optimizers */
  /*  might hold one of these in a register. */

  /*  Arguments */
  /*  ========= */

  /*  A, B    (input) DOUBLE PRECISION */
  /*          The values A and B. */

  /* ===================================================================== */

  /*     .. Executable Statements .. */

  ret_val = *a + *b;

  return ret_val;

  /*     End of DLAMC3 */

} /* dlamc3_ */


/* *********************************************************************** */

static int
dlamc4_(int * emin, double * start, int * base)
{
  /* System generated locals */
  int i__1;
  double d__1;

  /* Local variables */
  static double a;
  static int i__;
  static double b1, b2, c1, c2, d1, d2, one, zero, rbase;


  /*  -- LAPACK auxiliary routine (version 3.0) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLAMC4 is a service routine for DLAMC2. */

  /*  Arguments */
  /*  ========= */

  /*  EMIN    (output) EMIN */
  /*          The minimum exponent before (gradual) underflow, computed by */
  /*          setting A = START and dividing by BASE until the previous A */
  /*          can not be recovered. */

  /*  START   (input) DOUBLE PRECISION */
  /*          The starting point for determining EMIN. */

  /*  BASE    (input) INTEGER */
  /*          The base of the machine. */

  /* ===================================================================== */

  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  a = *start;
  one = 1.;
  rbase = one / *base;
  zero = 0.;
  *emin = 1;
  d__1 = a * rbase;
  b1 = dlamc3_(&d__1, &zero);
  c1 = a;
  c2 = a;
  d1 = a;
  d2 = a;
  /* +    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND. */
  /*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP */
L10:
  if (c1 == a && c2 == a && d1 == a && d2 == a)
  {
    --(*emin);
    a = b1;
    d__1 = a / *base;
    b1 = dlamc3_(&d__1, &zero);
    d__1 = b1 * *base;
    c1 = dlamc3_(&d__1, &zero);
    d1 = zero;
    i__1 = *base;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      d1 += b1;
      /* L20: */
    }
    d__1 = a * rbase;
    b2 = dlamc3_(&d__1, &zero);
    d__1 = b2 / rbase;
    c2 = dlamc3_(&d__1, &zero);
    d2 = zero;
    i__1 = *base;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
      d2 += b2;
      /* L30: */
    }
    goto L10;
  }
  /* +    END WHILE */

  return 0;

  /*     End of DLAMC4 */

} /* dlamc4_ */


/* *********************************************************************** */

static int
dlamc5_(int * beta, int * p, int * emin, bool * ieee, int * emax, double * rmax)
{
  /* System generated locals */
  int i__1;
  double d__1;

  /* Local variables */
  static int i__;
  static double y, z__;
  static int try__, lexp;
  static double oldy;
  static int uexp, nbits;
  static double recbas;
  static int exbits, expsum;


  /*  -- LAPACK auxiliary routine (version 3.0) -- */
  /*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
  /*     Courant Institute, Argonne National Lab, and Rice University */
  /*     October 31, 1992 */

  /*     .. Scalar Arguments .. */
  /*     .. */

  /*  Purpose */
  /*  ======= */

  /*  DLAMC5 attempts to compute RMAX, the largest machine floating-point */
  /*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum */
  /*  approximately to a power of 2.  It will fail on machines where this */
  /*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625, */
  /*  EMAX = 28718).  It will also fail if the value supplied for EMIN is */
  /*  too large (i.e. too close to zero), probably with overflow. */

  /*  Arguments */
  /*  ========= */

  /*  BETA    (input) INTEGER */
  /*          The base of floating-point arithmetic. */

  /*  P       (input) INTEGER */
  /*          The number of base BETA digits in the mantissa of a */
  /*          floating-point value. */

  /*  EMIN    (input) INTEGER */
  /*          The minimum exponent before (gradual) underflow. */

  /*  IEEE    (input) LOGICAL */
  /*          A bool flag specifying whether or not the arithmetic */
  /*          system is thought to comply with the IEEE standard. */

  /*  EMAX    (output) INTEGER */
  /*          The largest exponent before overflow */

  /*  RMAX    (output) DOUBLE PRECISION */
  /*          The largest machine floating-point number. */

  /* ===================================================================== */

  /*     .. Parameters .. */
  /*     .. */
  /*     .. Local Scalars .. */
  /*     .. */
  /*     .. External Functions .. */
  /*     .. */
  /*     .. Intrinsic Functions .. */
  /*     .. */
  /*     .. Executable Statements .. */

  /*     First compute LEXP and UEXP, two powers of 2 that bound */
  /*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum */
  /*     approximately to the bound that is closest to abs(EMIN). */
  /*     (EMAX is the exponent of the required number RMAX). */

  lexp = 1;
  exbits = 1;
L10:
  try__ = lexp << 1;
  if (try__ <= -(*emin))
  {
    lexp = try__;
    ++exbits;
    goto L10;
  }
  if (lexp == -(*emin))
  {
    uexp = lexp;
  }
  else
  {
    uexp = try__;
    ++exbits;
  }

  /*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater */
  /*     than or equal to EMIN. EXBITS is the number of bits needed to */
  /*     store the exponent. */

  if (uexp + *emin > -lexp - *emin)
  {
    expsum = lexp << 1;
  }
  else
  {
    expsum = uexp << 1;
  }

  /*     EXPSUM is the exponent range, approximately equal to */
  /*     EMAX - EMIN + 1 . */

  *emax = expsum + *emin - 1;
  nbits = exbits + 1 + *p;

  /*     NBITS is the total number of bits needed to store a */
  /*     floating-point number. */

  if (nbits % 2 == 1 && *beta == 2)
  {

    /*        Either there are an odd number of bits used to store a */
    /*        floating-point number, which is unlikely, or some bits are */
    /*        not used in the representation of numbers, which is possible, */
    /*        (e.g. Cray machines) or the mantissa has an implicit bit, */
    /*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the */
    /*        most likely. We have to assume the last alternative. */
    /*        If this is true, then we need to reduce EMAX by one because */
    /*        there must be some way of representing zero in an implicit-bit */
    /*        system. On machines like Cray, we are reducing EMAX by one */
    /*        unnecessarily. */

    --(*emax);
  }

  if (*ieee)
  {

    /*        Assume we are on an IEEE machine which reserves one exponent */
    /*        for infinity and NaN. */

    --(*emax);
  }

  /*     Now create RMAX, the largest machine number, which should */
  /*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX . */

  /*     First compute 1.0 - BETA**(-P), being careful that the */
  /*     result is less than 1.0 . */

  recbas = 1. / *beta;
  z__ = *beta - 1.;
  y = 0.;
  i__1 = *p;
  for (i__ = 1; i__ <= i__1; ++i__)
  {
    z__ *= recbas;
    if (y < 1.)
    {
      oldy = y;
    }
    y = dlamc3_(&y, &z__);
    /* L20: */
  }
  if (y >= 1.)
  {
    y = oldy;
  }

  /*     Now multiply by BETA**EMAX to get RMAX. */

  i__1 = *emax;
  for (i__ = 1; i__ <= i__1; ++i__)
  {
    d__1 = y * *beta;
    y = dlamc3_(&d__1, &c_b507);
    /* L30: */
  }

  *rmax = y;
  return 0;

  /*     End of DLAMC5 */

} /* dlamc5_ */
